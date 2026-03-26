# calcul_distances_adn.py

import os
import csv
import gc
import re
import time
import itertools
from datetime import datetime, date

# ============================================================
# Section 0 — Paramètres configurables
# ============================================================

CHUNK_SIZE     = 500
PAUSE_SECS     = 2
SEUIL_LIEN     = 0.01
TAUX_EVOLUTION = 2.2e-6   # substitutions/site/jour (calibré SARS-CoV-2)
BASES_AMBIGUES = set("N-?.")

COLONNES = [
    "id_seq1", "id_seq2",
    "pays_seq1", "pays_seq2",
    "trimestre_seq1", "trimestre_seq2",
    "hamming_brut", "longueur_ref", "nb_sites_valides",
    "dist_normalisee", "delta_jours", "dist_corrigee", "lien_possible"
]

# ============================================================
# Section 1 — Paramètres via input
# ============================================================

print("=" * 60)
print("  Calcul des distances évolutives ADN")
print("=" * 60)

FASTA_DIR    = input("\nDossier des fichiers FASTA (ex: data_merged_v4/sequences) : ").strip()
CANCER_DIR   = input("Dossier cancer   (ex: data_merged_v4/cancer)              : ").strip()
BIO_DIR      = input("Dossier biomedical (ex: data_merged_v4/biomedical)        : ").strip()
SOCIAL_DIR   = input("Dossier social   (ex: data_merged_v4/social)              : ").strip()

OUTPUT_FILE     = input("\nFichier de résultats (défaut: distances_output.csv)   : ").strip() or "distances_output.csv"
CHECKPOINT_FILE = input("Fichier checkpoint  (défaut: distances_checkpoint.csv) : ").strip() or "distances_checkpoint.csv"

CATEGORIES_DIRS = {
    "cancer":     CANCER_DIR,
    "biomedical": BIO_DIR,
    "social":     SOCIAL_DIR,
}

# Vérification que le dossier FASTA existe
if not os.path.isdir(FASTA_DIR):
    print(f"\n❌ Dossier FASTA introuvable : {FASTA_DIR}")
    exit(1)

# ============================================================
# Section 2 — Fonctions utilitaires
# ============================================================

def extraire_pays_trimestre(nom_fichier):
    """
    Parse un nom de fichier du type FR_2023Q1.fasta
    Retourne (pays, trimestre) ou (None, None).
    """
    base = os.path.splitext(nom_fichier)[0]
    m = re.match(r'^([A-Z]{2,3})_(\d{4}Q[1-4])$', base)
    if m:
        return m.group(1), m.group(2)
    return None, None


def parse_date_header(header):
    """
    Extrait la date depuis un header FASTA : >ID.jour.mois.année.pays
    Retourne un objet date ou None.
    """
    parts = header[1:].split(".")
    if len(parts) >= 4:
        try:
            return date(int(parts[3]), int(parts[2]), int(parts[1]))
        except (ValueError, IndexError):
            pass
    return None


def lire_fasta(path):
    """
    Lit un fichier FASTA.
    Retourne une liste de dicts { id, sequence, date_diagnostic }.
    """
    sequences = []
    header    = None
    seq       = ""

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequences.append({
                        "id":              header[1:].split(".")[0],
                        "sequence":        seq.upper(),
                        "date_diagnostic": parse_date_header(header),
                    })
                header = line
                seq    = ""
            else:
                seq += line

    if header is not None:
        sequences.append({
            "id":              header[1:].split(".")[0],
            "sequence":        seq.upper(),
            "date_diagnostic": parse_date_header(header),
        })

    return sequences

# ============================================================
# Section 3 — Chargement des données
# ============================================================

def charger_tous_fasta(fasta_dir):
    """
    Charge tous les .fasta du dossier.
    Retourne une liste de dicts { id, sequence, date_diagnostic, pays, trimestre }.
    """
    fichiers = sorted(f for f in os.listdir(fasta_dir) if f.endswith(".fasta"))
    print(f"\n  Chargement de {len(fichiers)} fichier(s) FASTA...")

    toutes = []
    for fname in fichiers:
        pays, trimestre = extraire_pays_trimestre(fname)
        if pays is None:
            print(f"  ⚠  Nom ignoré (format inattendu) : {fname}")
            continue

        seqs = lire_fasta(os.path.join(fasta_dir, fname))
        for s in seqs:
            s["pays"]      = pays
            s["trimestre"] = trimestre
        toutes.extend(seqs)

    print(f"  → {len(toutes)} séquences chargées")
    return toutes


def charger_tous_metadata(categories_dirs):
    """
    Charge tous les CSV depuis les dossiers catégorie/pays/*.csv.
    Retourne un dict { id → date_la_plus_precoce }.
    """
    dates_par_id = {}

    for categorie, cat_dir in categories_dirs.items():
        if not os.path.isdir(cat_dir):
            print(f"  ⚠  Dossier absent, ignoré : {cat_dir}")
            continue

        for pays in os.listdir(cat_dir):
            pays_dir = os.path.join(cat_dir, pays)
            if not os.path.isdir(pays_dir):
                continue

            for fname in os.listdir(pays_dir):
                if not fname.endswith(".csv"):
                    continue

                path = os.path.join(pays_dir, fname)
                with open(path, "r", encoding="utf-8", newline="") as f:
                    reader = csv.DictReader(f)
                    fields = reader.fieldnames or []

                    if "id" not in fields:
                        print(f"  ⚠  Colonne 'id' absente : {path}")
                        continue

                    has_date = "date" in fields
                    if not has_date:
                        print(f"  ⚠  Colonne 'date' absente : {path}")

                    for row in reader:
                        row_id = row.get("id", "").strip()
                        if not row_id:
                            continue

                        date_val = None
                        if has_date and row.get("date", "").strip():
                            try:
                                date_val = datetime.strptime(
                                    row["date"].strip(), "%Y-%m-%d"
                                ).date()
                            except ValueError:
                                pass

                        # Garde la date la plus précoce si ID déjà vu
                        if row_id not in dates_par_id:
                            dates_par_id[row_id] = date_val
                        elif date_val is not None:
                            if dates_par_id[row_id] is None or date_val < dates_par_id[row_id]:
                                dates_par_id[row_id] = date_val

    print(f"  → {len(dates_par_id)} ID avec métadonnées")
    return dates_par_id


def fusionner_dates(sequences, dates_meta):
    """
    Complète la date manquante d'une séquence avec les métadonnées CSV.
    Priorité : date du header FASTA, sinon date du CSV.
    """
    for s in sequences:
        if s["date_diagnostic"] is None:
            s["date_diagnostic"] = dates_meta.get(s["id"])
    return sequences

# ============================================================
# Section 4 — Calcul des distances
# ============================================================

def hamming_brut(seq1, seq2):
    """
    Calcule la distance de Hamming en ignorant les bases ambiguës.
    Retourne (hamming, longueur_max, nb_sites_valides).
    """
    longueur = max(len(seq1), len(seq2))
    hamming  = 0
    valides  = 0

    for i in range(min(len(seq1), len(seq2))):
        b1, b2 = seq1[i], seq2[i]
        if b1 in BASES_AMBIGUES or b2 in BASES_AMBIGUES:
            continue
        valides += 1
        if b1 != b2:
            hamming += 1

    return hamming, longueur, valides


def corriger_distance(dist_norm, date1, date2):
    """
    Applique la correction temporelle :
    dist_corrigee = max(0, dist_norm - TAUX_EVOLUTION × delta_jours)
    Retourne (dist_corrigee, delta_jours) — None si date manquante.
    """
    if dist_norm is None or date1 is None or date2 is None:
        return None, None
    delta_jours = abs((date2 - date1).days)
    dist_corr   = max(0.0, dist_norm - TAUX_EVOLUTION * delta_jours)
    return dist_corr, delta_jours

# ============================================================
# Section 5 — Checkpoint
# ============================================================

def charger_checkpoint(checkpoint_file):
    """Retourne un set de clés 'id1|||id2' déjà calculées."""
    deja_faites = set()
    if not os.path.isfile(checkpoint_file):
        return deja_faites

    with open(checkpoint_file, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            deja_faites.add(f"{row['id_seq1']}|||{row['id_seq2']}")

    print(f"  ♻  Reprise : {len(deja_faites):,} paire(s) déjà calculée(s)")
    return deja_faites


def init_checkpoint(checkpoint_file):
    """Crée le checkpoint avec l'en-tête si absent."""
    if not os.path.isfile(checkpoint_file):
        with open(checkpoint_file, "w", encoding="utf-8", newline="") as f:
            csv.DictWriter(f, fieldnames=COLONNES).writeheader()


def append_checkpoint(checkpoint_file, lignes):
    """Ajoute des lignes au checkpoint sans réécrire l'en-tête."""
    with open(checkpoint_file, "a", encoding="utf-8", newline="") as f:
        csv.DictWriter(f, fieldnames=COLONNES).writerows(lignes)

# ============================================================
# Section 6 — Calcul par chunks
# ============================================================

def calculer_paires(sequences, deja_faites, checkpoint_file):
    n              = len(sequences)
    nb_total       = n * (n - 1) // 2
    nb_restantes   = nb_total - len(deja_faites)
    nb_chunks      = max(1, (nb_restantes + CHUNK_SIZE - 1) // CHUNK_SIZE)

    print(f"\n  {n} séquences → {nb_total:,} paires totales")
    print(f"  {nb_restantes:,} restante(s) → {nb_chunks} chunk(s) de {CHUNK_SIZE}")

    chunk        = []
    nb_calculees = 0
    chunk_num    = 0

    for i, j in itertools.combinations(range(n), 2):
        s1  = sequences[i]
        s2  = sequences[j]
        cle = f"{s1['id']}|||{s2['id']}"

        if cle in deja_faites:
            continue

        h, lon, val   = hamming_brut(s1["sequence"], s2["sequence"])
        dist_norm     = (h / val) if val > 0 else None
        dist_corr, dj = corriger_distance(dist_norm, s1["date_diagnostic"], s2["date_diagnostic"])
        lien          = (dist_corr <= SEUIL_LIEN) if dist_corr is not None else None

        chunk.append({
            "id_seq1":          s1["id"],
            "id_seq2":          s2["id"],
            "pays_seq1":        s1["pays"],
            "pays_seq2":        s2["pays"],
            "trimestre_seq1":   s1["trimestre"],
            "trimestre_seq2":   s2["trimestre"],
            "hamming_brut":     h,
            "longueur_ref":     lon,
            "nb_sites_valides": val,
            "dist_normalisee":  round(dist_norm, 8) if dist_norm is not None else "",
            "delta_jours":      dj if dj is not None else "",
            "dist_corrigee":    round(dist_corr, 8) if dist_corr is not None else "",
            "lien_possible":    lien if lien is not None else "",
        })

        if len(chunk) >= CHUNK_SIZE:
            chunk_num    += 1
            nb_calculees += len(chunk)
            append_checkpoint(checkpoint_file, chunk)
            pct = nb_calculees / nb_restantes * 100 if nb_restantes else 100
            print(f"  Chunk {chunk_num}/{nb_chunks} — "
                  f"{nb_calculees:,}/{nb_restantes:,} ({pct:.1f}%)")
            chunk = []
            gc.collect()
            time.sleep(PAUSE_SECS)

    # Dernier chunk
    if chunk:
        chunk_num    += 1
        nb_calculees += len(chunk)
        append_checkpoint(checkpoint_file, chunk)
        print(f"  Chunk {chunk_num}/{nb_chunks} — "
              f"{nb_calculees:,}/{nb_restantes:,} (100%)")
        gc.collect()

    return nb_calculees

# ============================================================
# Section 7 — Export final
# ============================================================

def export_final(checkpoint_file, output_file):
    lignes = []
    with open(checkpoint_file, "r", encoding="utf-8", newline="") as f:
        lignes = list(csv.DictReader(f))

    lignes.sort(key=lambda r: (r["id_seq1"], r["id_seq2"]))

    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLONNES)
        writer.writeheader()
        writer.writerows(lignes)

    nb_liens = sum(
        1 for r in lignes
        if str(r.get("lien_possible", "")).strip().lower() == "true"
    )

    print(f"\n{'='*60}")
    print(f"  ✅ Terminé")
    print(f"  Total paires           : {len(lignes):,}")
    print(f"  Liens possibles        : {nb_liens:,}")
    print(f"  Fichier de résultats   : {output_file}")
    print(f"{'='*60}")

# ============================================================
# Main
# ============================================================

print("\n" + "=" * 60)
print("  Résumé")
print("=" * 60)
print(f"  FASTA      : {FASTA_DIR}")
for cat, d in CATEGORIES_DIRS.items():
    print(f"  {cat:12}: {d}")
print(f"  Checkpoint : {CHECKPOINT_FILE}")
print(f"  Sortie     : {OUTPUT_FILE}")
print(f"  Seuil lien : {SEUIL_LIEN}  |  Taux évolution : {TAUX_EVOLUTION} sub/site/jour")
print("=" * 60)

confirm = input("\nLancer le calcul ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)

# — Chargement —
print("\n─" * 30)
print("Chargement FASTA")
sequences = charger_tous_fasta(FASTA_DIR)

print("\nChargement métadonnées CSV")
dates_meta = charger_tous_metadata(CATEGORIES_DIRS)

sequences = fusionner_dates(sequences, dates_meta)
sans_date = sum(1 for s in sequences if s["date_diagnostic"] is None)
if sans_date:
    print(f"  ⚠  {sans_date} séquence(s) sans date — dist_corrigee sera vide pour ces paires")

# — Checkpoint —
print("\nVérification checkpoint")
init_checkpoint(CHECKPOINT_FILE)
deja_faites = charger_checkpoint(CHECKPOINT_FILE)

# — Calcul —
print("\nCalcul des distances")
debut        = time.time()
nb_calculees = calculer_paires(sequences, deja_faites, CHECKPOINT_FILE)
duree        = time.time() - debut
print(f"\n  Temps : {duree:.1f}s ({duree/60:.1f} min)")

# — Export —
print("\nExport final")
export_final(CHECKPOINT_FILE, OUTPUT_FILE)

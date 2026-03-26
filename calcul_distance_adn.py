# calcul_distances_adn.py

import os
import csv
import gc
import time
import itertools
import re
from datetime import datetime, date

# ============================================================
# Section 0 — Paramètres configurables
# ============================================================

RACINE          = "data_merged_v4"
OUTPUT_FILE     = "distances_output.csv"
CHECKPOINT_FILE = "distances_checkpoint.csv"
CHUNK_SIZE      = 500
PAUSE_SECS      = 2
SEUIL_LIEN      = 0.01
TAUX_EVOLUTION  = 2.2e-6   # substitutions/site/jour (calibré SARS-CoV-2)
CATEGORIES_META = ["cancer", "biomedical", "social"]

BASES_AMBIGUES  = set("N-?.")

# Colonnes du fichier de sortie
COLONNES = [
    "id_seq1", "id_seq2",
    "pays_seq1", "pays_seq2",
    "trimestre_seq1", "trimestre_seq2",
    "hamming_brut", "longueur_ref", "nb_sites_valides",
    "dist_normalisee", "delta_jours", "dist_corrigee", "lien_possible"
]

# ============================================================
# Section 1 — Fonctions utilitaires
# ============================================================

def extraire_pays_trimestre(nom_fichier):
    """
    Parse un nom de fichier FASTA/CSV du type CP_XXXXQX.
    Retourne (pays, trimestre) ou (None, None) si le format ne correspond pas.
    Exemple : FR_2023Q1.fasta → ("FR", "2023Q1")
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
            jour  = int(parts[1])
            mois  = int(parts[2])
            annee = int(parts[3])
            return date(annee, mois, jour)
        except (ValueError, IndexError):
            pass
    return None


def lire_fasta(path):
    """
    Lit un fichier FASTA.
    Retourne une liste de dicts :
      { id, sequence, date_diagnostic }
    La date est extraite directement depuis le header.
    """
    sequences = []
    header = None
    seq    = ""

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq_id = header[1:].split(".")[0]
                    sequences.append({
                        "id":               seq_id,
                        "sequence":         seq.upper(),
                        "date_diagnostic":  parse_date_header(header),
                    })
                header = line
                seq    = ""
            else:
                seq += line

    if header is not None:
        seq_id = header[1:].split(".")[0]
        sequences.append({
            "id":              seq_id,
            "sequence":        seq.upper(),
            "date_diagnostic": parse_date_header(header),
        })

    return sequences


# ============================================================
# Section 2 — Chargement des données
# ============================================================

def charger_tous_fasta(racine):
    """
    Scanne racine/sequences/, charge tous les .fasta.
    Retourne une liste de dicts { id, sequence, date_diagnostic, pays, trimestre }.
    """
    seq_dir = os.path.join(racine, "sequences")
    if not os.path.isdir(seq_dir):
        raise FileNotFoundError(f"Dossier introuvable : {seq_dir}")

    toutes = []
    fichiers = sorted(f for f in os.listdir(seq_dir) if f.endswith(".fasta"))

    print(f"  Chargement de {len(fichiers)} fichier(s) FASTA...")

    for fname in fichiers:
        pays, trimestre = extraire_pays_trimestre(fname)
        if pays is None:
            print(f"  ⚠  Nom de fichier ignoré (format inattendu) : {fname}")
            continue

        seqs = lire_fasta(os.path.join(seq_dir, fname))
        for s in seqs:
            s["pays"]      = pays
            s["trimestre"] = trimestre
        toutes.extend(seqs)

    print(f"  → {len(toutes)} séquences chargées")
    return toutes


def charger_tous_metadata(racine, categories):
    """
    Scanne racine/{categorie}/{pays}/*.csv, charge tous les CSV.
    Retourne un dict { id → date } en gardant la date la plus précoce
    si un ID apparaît dans plusieurs catégories.
    Attend les colonnes 'id' et 'date'.
    """
    dates_par_id = {}  # id → date la plus précoce

    for categorie in categories:
        cat_dir = os.path.join(racine, categorie)
        if not os.path.isdir(cat_dir):
            print(f"  ⚠  Dossier absent : {cat_dir}")
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
                    if "id" not in (reader.fieldnames or []):
                        print(f"  ⚠  Colonne 'id' absente : {path}")
                        continue
                    has_date = "date" in (reader.fieldnames or [])
                    if not has_date:
                        print(f"  ⚠  Colonne 'date' absente : {path} — dates NA")

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

                        # On garde la date la plus précoce
                        if row_id not in dates_par_id:
                            dates_par_id[row_id] = date_val
                        elif date_val is not None:
                            if dates_par_id[row_id] is None or date_val < dates_par_id[row_id]:
                                dates_par_id[row_id] = date_val

    print(f"  → {len(dates_par_id)} ID avec métadonnées chargés")
    return dates_par_id


# ============================================================
# Section 3 — Fonctions de calcul
# ============================================================

def hamming_brut(seq1, seq2):
    """
    Calcule la distance de Hamming entre deux séquences.
    Ignore les positions ambiguës (N, -, ?, .).
    Retourne (hamming, longueur_max, nb_sites_valides).
    """
    longueur = max(len(seq1), len(seq2))
    hamming  = 0
    valides  = 0

    for i in range(min(len(seq1), len(seq2))):
        b1 = seq1[i]
        b2 = seq2[i]
        if b1 in BASES_AMBIGUES or b2 in BASES_AMBIGUES:
            continue
        valides += 1
        if b1 != b2:
            hamming += 1

    return hamming, longueur, valides


def hamming_normalise(hamming, nb_valides):
    """Retourne hamming / nb_valides, ou None si nb_valides == 0."""
    if nb_valides == 0:
        return None
    return hamming / nb_valides


def corriger_distance(dist_norm, date1, date2):
    """
    Soustrait la dérive évolutive attendue entre deux dates.
    dist_corrigee = max(0, dist_norm - TAUX_EVOLUTION × delta_jours)
    Retourne (dist_corrigee, delta_jours) — les deux peuvent être None
    si l'une des dates est manquante.
    """
    if dist_norm is None:
        return None, None
    if date1 is None or date2 is None:
        return None, None

    delta_jours = abs((date2 - date1).days)
    derive      = TAUX_EVOLUTION * delta_jours
    dist_corr   = max(0.0, dist_norm - derive)
    return dist_corr, delta_jours


# ============================================================
# Section 4 — Fusion séquences + dates
# ============================================================

def fusionner_dates(sequences, dates_meta):
    """
    Pour chaque séquence, cherche sa date dans les métadonnées.
    Priorité : date du header FASTA, sinon date des métadonnées CSV.
    """
    for s in sequences:
        if s["date_diagnostic"] is None:
            s["date_diagnostic"] = dates_meta.get(s["id"])
    return sequences


# ============================================================
# Section 5 — Reprise sur crash
# ============================================================

def charger_checkpoint(checkpoint_file):
    """
    Lit les paires déjà calculées depuis le checkpoint.
    Retourne un set de clés "id1|||id2".
    """
    deja_faites = set()
    if not os.path.isfile(checkpoint_file):
        return deja_faites

    with open(checkpoint_file, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            cle = f"{row['id_seq1']}|||{row['id_seq2']}"
            deja_faites.add(cle)

    print(f"  ♻  Reprise : {len(deja_faites)} paire(s) déjà calculée(s)")
    return deja_faites


def init_checkpoint(checkpoint_file):
    """Crée le fichier checkpoint avec l'en-tête si absent."""
    if not os.path.isfile(checkpoint_file):
        with open(checkpoint_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=COLONNES)
            writer.writeheader()


def append_checkpoint(checkpoint_file, lignes):
    """Ajoute des lignes au checkpoint sans réécrire l'en-tête."""
    with open(checkpoint_file, "a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLONNES)
        writer.writerows(lignes)


# ============================================================
# Section 6 — Calcul par chunks
# ============================================================

def calculer_paires(sequences, deja_faites, checkpoint_file):
    """
    Génère toutes les paires (triangulaire supérieure), filtre celles
    déjà calculées, puis calcule par chunks de CHUNK_SIZE.
    """
    n = len(sequences)
    nb_paires_total = n * (n - 1) // 2
    print(f"\n  {n} séquences → {nb_paires_total:,} paires totales")

    # Génère les paires non encore calculées
    def paires_restantes():
        for i, j in itertools.combinations(range(n), 2):
            s1 = sequences[i]
            s2 = sequences[j]
            cle = f"{s1['id']}|||{s2['id']}"
            if cle not in deja_faites:
                yield s1, s2

    nb_restantes  = nb_paires_total - len(deja_faites)
    nb_chunks     = (nb_restantes + CHUNK_SIZE - 1) // CHUNK_SIZE
    print(f"  {nb_restantes:,} paire(s) restante(s) → {nb_chunks} chunk(s) de {CHUNK_SIZE}")

    chunk       = []
    nb_calculees = 0
    chunk_num   = 0

    for s1, s2 in paires_restantes():
        h, lon, val     = hamming_brut(s1["sequence"], s2["sequence"])
        dist_norm       = hamming_normalise(h, val)
        dist_corr, dj   = corriger_distance(
            dist_norm, s1["date_diagnostic"], s2["date_diagnostic"]
        )

        lien = None
        if dist_corr is not None:
            lien = dist_corr <= SEUIL_LIEN

        chunk.append({
            "id_seq1":         s1["id"],
            "id_seq2":         s2["id"],
            "pays_seq1":       s1["pays"],
            "pays_seq2":       s2["pays"],
            "trimestre_seq1":  s1["trimestre"],
            "trimestre_seq2":  s2["trimestre"],
            "hamming_brut":    h,
            "longueur_ref":    lon,
            "nb_sites_valides":val,
            "dist_normalisee": round(dist_norm, 8) if dist_norm is not None else "",
            "delta_jours":     dj if dj is not None else "",
            "dist_corrigee":   round(dist_corr, 8) if dist_corr is not None else "",
            "lien_possible":   lien if lien is not None else "",
        })

        if len(chunk) >= CHUNK_SIZE:
            chunk_num   += 1
            nb_calculees += len(chunk)
            append_checkpoint(checkpoint_file, chunk)
            pct = nb_calculees / nb_restantes * 100 if nb_restantes else 100
            print(f"  Chunk {chunk_num}/{nb_chunks} — "
                  f"{nb_calculees:,}/{nb_restantes:,} paires ({pct:.1f}%)")
            chunk = []
            gc.collect()
            time.sleep(PAUSE_SECS)

    # Dernier chunk incomplet
    if chunk:
        chunk_num   += 1
        nb_calculees += len(chunk)
        append_checkpoint(checkpoint_file, chunk)
        print(f"  Chunk {chunk_num}/{nb_chunks} — "
              f"{nb_calculees:,}/{nb_restantes:,} paires (100%)")
        gc.collect()

    return nb_calculees


# ============================================================
# Section 7 — Export final
# ============================================================

def export_final(checkpoint_file, output_file):
    """
    Lit le checkpoint complet, trie par id_seq1 puis id_seq2,
    exporte dans output_file, affiche un résumé.
    """
    lignes = []
    with open(checkpoint_file, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lignes.append(row)

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
    print(f"  ✅ Export final terminé")
    print(f"  Total paires  : {len(lignes):,}")
    print(f"  Liens possibles détectés : {nb_liens:,}")
    print(f"  Fichier de résultats : {output_file}")
    print(f"{'='*60}")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("  Calcul des distances évolutives ADN")
    print("=" * 60)
    print(f"\n  Dossier racine   : {RACINE}")
    print(f"  Seuil lien       : {SEUIL_LIEN}")
    print(f"  Taux évolution   : {TAUX_EVOLUTION} sub/site/jour")
    print(f"  Chunk size       : {CHUNK_SIZE}")
    print(f"  Checkpoint       : {CHECKPOINT_FILE}")
    print()

    # — Chargement —
    print("─" * 60)
    print("  Chargement des séquences FASTA")
    print("─" * 60)
    sequences = charger_tous_fasta(RACINE)

    print()
    print("─" * 60)
    print("  Chargement des métadonnées CSV")
    print("─" * 60)
    dates_meta = charger_tous_metadata(RACINE, CATEGORIES_META)

    # — Fusion dates —
    sequences = fusionner_dates(sequences, dates_meta)
    sans_date = sum(1 for s in sequences if s["date_diagnostic"] is None)
    if sans_date:
        print(f"\n  ⚠  {sans_date} séquence(s) sans date — "
              f"dist_corrigee sera NA pour ces paires")

    # — Reprise —
    print()
    print("─" * 60)
    print("  Vérification du checkpoint")
    print("─" * 60)
    init_checkpoint(CHECKPOINT_FILE)
    deja_faites = charger_checkpoint(CHECKPOINT_FILE)

    # — Calcul —
    print()
    print("─" * 60)
    print("  Calcul des distances")
    print("─" * 60)
    debut = time.time()
    nb_calculees = calculer_paires(sequences, deja_faites, CHECKPOINT_FILE)
    duree = time.time() - debut
    print(f"\n  Temps de calcul : {duree:.1f}s "
          f"({duree/60:.1f} min)")

    # — Export —
    print()
    print("─" * 60)
    print("  Export final")
    print("─" * 60)
    export_final(CHECKPOINT_FILE, OUTPUT_FILE)

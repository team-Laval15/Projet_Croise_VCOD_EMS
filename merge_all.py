# merge_all.py

import os
import csv
import io
from pipelines.utils import ensure_dir

# ============================================================
# Paramètres
# ============================================================

print("=" * 60)
print("     Fusion des données — détection de doublons")
print("=" * 60)

FASTA_OLD_DIR = input("\nDossier anciens FASTA (ex: fasta_data)     : ").strip() or "fasta_data"
FASTA_NEW_DIR = input("Dossier nouveaux FASTA (ex: sequences)     : ").strip() or "sequences"

CSV_OLD_DIR   = input("\nDossier anciens CSV (ex: csv_downloads)    : ").strip() or "csv_downloads"
CSV_NEW_DIR   = input("Dossier nouveaux CSV (ex: data_final)      : ").strip() or "data_final"

OUT_DIR       = input("\nDossier de sortie fusionné (ex: merged)    : ").strip() or "merged"
DUP_DIR       = input("Dossier pour les doublons (ex: duplicates) : ").strip() or "duplicates"

CSV_DATASETS  = ["cancer", "biomedical", "social"]

# Création des dossiers de sortie
ensure_dir(os.path.join(OUT_DIR, "sequences"))
ensure_dir(os.path.join(DUP_DIR, "sequences"))

# ============================================================
# Résumé
# ============================================================

print("\n" + "=" * 60)
print("  Résumé")
print("=" * 60)
print(f"  FASTA anciens  : {FASTA_OLD_DIR}/")
print(f"  FASTA nouveaux : {FASTA_NEW_DIR}/")
print(f"  CSV anciens    : {CSV_OLD_DIR}/{{dataset}}/{{pays}}/")
print(f"  CSV nouveaux   : {CSV_NEW_DIR}/{{dataset}}/{{pays}}/")
print(f"  Sortie fusion  : {OUT_DIR}/")
print(f"  Doublons       : {DUP_DIR}/")
print(f"  Structure CSV  : {{dataset}}/{{pays}}/{{fichier}}.csv")
print("=" * 60)

confirm = input("\nLancer la fusion ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)

print()

# ============================================================
# Helpers généraux
# ============================================================

def list_files(folder, extension):
    """Liste tous les fichiers avec l'extension donnée dans un dossier."""
    if not os.path.isdir(folder):
        return []
    return [f for f in os.listdir(folder) if f.endswith(extension)]


def read_fasta(path):
    """Lit un fichier FASTA, retourne une liste de (header, sequence)."""
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
                    sequences.append((header, seq))
                header = line
                seq    = ""
            else:
                seq += line
    if header is not None:
        sequences.append((header, seq))
    return sequences


def write_fasta(path, sequences):
    """Écrit une liste de (header, seq) dans un fichier FASTA."""
    with open(path, "w", encoding="utf-8") as f:
        for header, seq in sequences:
            f.write(f"{header}\n{seq}\n")


def extract_fasta_id(header):
    """Extrait l'ID depuis un header FASTA : >U58801.23.7.1987.FR → U58801"""
    return header[1:].split(".")[0]


def header_suffix(header):
    """Retourne tout ce qui suit le premier champ : .23.7.1987.FR"""
    parts = header[1:].split(".", 1)
    return "." + parts[1] if len(parts) > 1 else ""


# ============================================================
# Helpers génération de nouvel ID
# ============================================================

import re

def detect_id_format(existing_id):
    """
    Analyse un ID existant et retourne (prefix, digits, suffix) ou None.
    Exemples :
      U58801   → ("U",  "58801",  "")
      AF461909 → ("AF", "461909", "")
      AY010335 → ("AY", "010335", "")
      JQ292417 → ("JQ", "292417", "")
    """
    m = re.match(r'^([A-Za-z]+)(\d+)([A-Za-z]*)$', existing_id)
    if m:
        return m.group(1), m.group(2), m.group(3)
    return None


def generate_new_id(existing_ids):
    """
    Génère un nouvel ID unique en se basant sur le format des IDs existants.
    Prend l'ID le plus long/représentatif, incrémente le numéro jusqu'à
    trouver un ID absent de existing_ids.
    """
    # Choisit un ID de référence parmi ceux qui ont un format détectable
    ref_id = None
    ref_fmt = None
    for eid in existing_ids:
        fmt = detect_id_format(eid)
        if fmt:
            ref_id  = eid
            ref_fmt = fmt
            break

    if ref_fmt is None:
        # Fallback : format générique si aucun ID analysable
        n = 1
        while f"NEW{n:06d}" in existing_ids:
            n += 1
        return f"NEW{n:06d}"

    prefix, digits, suffix = ref_fmt
    num_len = len(digits)
    base    = int(digits)

    # Incrémente jusqu'à trouver un ID libre
    candidate = base + 1
    while True:
        new_id = f"{prefix}{str(candidate).zfill(num_len)}{suffix}"
        if new_id not in existing_ids:
            return new_id
        candidate += 1


# ============================================================
# Fusion FASTA
# ============================================================

print("─" * 60)
print("  Fusion FASTA")
print("─" * 60)

old_fasta_files = set(list_files(FASTA_OLD_DIR, ".fasta"))
new_fasta_files = set(list_files(FASTA_NEW_DIR, ".fasta"))
dup_fasta_files = set(list_files(os.path.join(DUP_DIR, "sequences"), ".fasta"))
all_fasta_files = old_fasta_files | new_fasta_files | dup_fasta_files

fasta_stats = {"fichiers": 0, "sequences": 0, "doublons_mis_de_cote": 0,
               "doublons_identiques_supprimes": 0, "doublons_reintegres": 0}

for filename in sorted(all_fasta_files):
    old_path = os.path.join(FASTA_OLD_DIR, filename)
    new_path = os.path.join(FASTA_NEW_DIR, filename)
    dup_path = os.path.join(DUP_DIR, "sequences", filename)
    out_path = os.path.join(OUT_DIR, "sequences", filename)

    # ── Étape 1 : lecture de toutes les sources ──────────────
    # Ordre de priorité : anciens → nouveaux
    sequences_vues = {}   # id → (header, seq) — première occurrence gardée
    candidats_dup  = []   # [(id_original, header, seq)] — doublons à analyser

    sources = []
    if filename in old_fasta_files:
        sources.append(old_path)
    if filename in new_fasta_files:
        sources.append(new_path)

    for src in sources:
        for header, seq in read_fasta(src):
            seq_id = extract_fasta_id(header)
            if seq_id not in sequences_vues:
                sequences_vues[seq_id] = (header, seq)
            else:
                candidats_dup.append((seq_id, header, seq))

    # ── Étape 2 : réintégration des doublons du run précédent ─
    # Les séquences dans duplicates/ ont déjà un nouvel ID → on les
    # réintègre directement si leur ID ne crée pas de nouveau conflit.
    reintegres = 0
    if filename in dup_fasta_files:
        for header, seq in read_fasta(dup_path):
            seq_id = extract_fasta_id(header)
            if seq_id not in sequences_vues:
                sequences_vues[seq_id] = (header, seq)
                reintegres += 1
            else:
                # Conflit encore présent → on le remet dans candidats
                candidats_dup.append((seq_id, header, seq))

    fasta_stats["doublons_reintegres"] += reintegres

    # ── Étape 3 : traitement des doublons restants ────────────
    # Pour chaque doublon, on compare la séquence avec l'originale :
    #   - identique → supprimé silencieusement
    #   - différente → nouvel ID généré et mis de côté
    nouveaux_doublons = []   # (new_header, seq) à écrire dans duplicates/
    supprimes         = 0

    for seq_id_original, header, seq in candidats_dup:
        _, seq_originale = sequences_vues[seq_id_original]

        if seq == seq_originale:
            # Séquences identiques → doublon réel, on supprime
            supprimes += 1
        else:
            # Séquences différentes → on génère un nouvel ID unique
            all_known_ids = set(sequences_vues.keys()) | {
                extract_fasta_id(h) for h, _ in nouveaux_doublons
            }
            new_id     = generate_new_id(all_known_ids)
            new_header = f">{new_id}{header_suffix(header)}"
            nouveaux_doublons.append((new_header, seq))

    fasta_stats["doublons_identiques_supprimes"] += supprimes
    fasta_stats["doublons_mis_de_cote"]          += len(nouveaux_doublons)

    # ── Étape 4 : écriture ───────────────────────────────────
    write_fasta(out_path, list(sequences_vues.values()))

    # Écrase l'ancien fichier duplicates avec les nouveaux doublons
    # (si vide ou inexistant on supprime le fichier pour ne pas laisser
    #  de résidu d'un run précédent)
    if nouveaux_doublons:
        write_fasta(dup_path, nouveaux_doublons)
    elif os.path.isfile(dup_path):
        os.remove(dup_path)

    fasta_stats["fichiers"]  += 1
    fasta_stats["sequences"] += len(sequences_vues)

    # ── Affichage ─────────────────────────────────────────────
    msg = f"  ✔  {filename}  ({len(sequences_vues)} séquences"
    details = []
    if supprimes:
        details.append(f"{supprimes} identique(s) supprimé(s)")
    if len(nouveaux_doublons):
        details.append(f"{len(nouveaux_doublons)} différent(s) mis de côté avec nouvel ID")
    if reintegres:
        details.append(f"{reintegres} doublon(s) réintégré(s)")
    if details:
        msg += " — " + ", ".join(details)
    print(msg + ")")

print(f"\n  → {fasta_stats['fichiers']} fichier(s) fusionné(s)")
print(f"     {fasta_stats['sequences']} séquences uniques dans merged/")
print(f"     {fasta_stats['doublons_identiques_supprimes']} doublon(s) identiques supprimés")
print(f"     {fasta_stats['doublons_mis_de_cote']} doublon(s) différents mis de côté (nouvel ID)")
if fasta_stats['doublons_reintegres']:
    print(f"     {fasta_stats['doublons_reintegres']} doublon(s) réintégrés depuis duplicates/")

# ============================================================
# Fusion CSV
# ============================================================

print()
print("─" * 60)
print("  Fusion CSV")
print("─" * 60)

csv_stats_total = {"fichiers": 0, "lignes": 0, "doublons": 0}

for dataset in CSV_DATASETS:
    old_ds_dir = os.path.join(CSV_OLD_DIR, dataset)
    new_ds_dir = os.path.join(CSV_NEW_DIR, dataset)

    # Collecte tous les pays présents dans l'un ou l'autre dossier
    old_countries = set(os.listdir(old_ds_dir)) if os.path.isdir(old_ds_dir) else set()
    new_countries = set(os.listdir(new_ds_dir)) if os.path.isdir(new_ds_dir) else set()
    all_countries = old_countries | new_countries

    if not all_countries:
        print(f"  ⚠  Aucun dossier pays trouvé pour [{dataset}]")
        continue

    ds_stats = {"fichiers": 0, "lignes": 0, "doublons": 0}
    print(f"\n  [{dataset}]")

    for country in sorted(all_countries):
        old_country_dir = os.path.join(old_ds_dir, country)
        new_country_dir = os.path.join(new_ds_dir, country)

        old_csv_files = set(list_files(old_country_dir, ".csv"))
        new_csv_files = set(list_files(new_country_dir, ".csv"))
        all_csv_files = old_csv_files | new_csv_files

        if not all_csv_files:
            continue

        # Création des dossiers de sortie pour ce pays
        ensure_dir(os.path.join(OUT_DIR, dataset, country))
        ensure_dir(os.path.join(DUP_DIR, dataset, country))

        print(f"    [{country}]")

        for filename in sorted(all_csv_files):
            old_path = os.path.join(old_country_dir, filename)
            new_path = os.path.join(new_country_dir, filename)

            rows_vues    = {}  # id → row dict — première occurrence gardée
            doublons_csv = {}  # id → row dict — occurrences en doublon
            fieldnames   = None

            # Lecture dans l'ordre : anciens d'abord, nouveaux ensuite
            sources = []
            if filename in old_csv_files:
                sources.append(old_path)
            if filename in new_csv_files:
                sources.append(new_path)

            for src in sources:
                with open(src, "r", encoding="utf-8", newline="") as f:
                    reader = csv.DictReader(f)
                    if fieldnames is None:
                        fieldnames = reader.fieldnames
                    for row in reader:
                        row_id = row.get("id")
                        if row_id is None:
                            fake_id = f"__noid_{len(rows_vues)}"
                            rows_vues[fake_id] = row
                        elif row_id not in rows_vues:
                            rows_vues[row_id] = row
                        else:
                            doublons_csv[row_id] = row

            if not fieldnames:
                print(f"      ⚠  {filename} — impossible de lire les colonnes, ignoré")
                continue

            # Écriture du fichier fusionné
            out_path = os.path.join(OUT_DIR, dataset, country, filename)
            with open(out_path, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(rows_vues.values())

            # Écriture des doublons
            if doublons_csv:
                dup_path = os.path.join(DUP_DIR, dataset, country, filename)
                with open(dup_path, "w", encoding="utf-8", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(doublons_csv.values())

            ds_stats["fichiers"]  += 1
            ds_stats["lignes"]    += len(rows_vues)
            ds_stats["doublons"]  += len(doublons_csv)

            status = f"      ✔  {filename}  ({len(rows_vues)} lignes"
            if doublons_csv:
                status += f", {len(doublons_csv)} doublon(s) mis de côté"
            print(status + ")")

    print(f"    → {ds_stats['fichiers']} fichier(s), "
          f"{ds_stats['lignes']} lignes uniques, "
          f"{ds_stats['doublons']} doublon(s)")

    csv_stats_total["fichiers"] += ds_stats["fichiers"]
    csv_stats_total["lignes"]   += ds_stats["lignes"]
    csv_stats_total["doublons"] += ds_stats["doublons"]

# ============================================================
# Fin
# ============================================================

print()
print("=" * 60)
print("  ✅ Fusion terminée")
print("=" * 60)
print(f"  FASTA : {fasta_stats['fichiers']} fichiers, "
      f"{fasta_stats['sequences']} séquences uniques, "
      f"{fasta_stats['doublons_identiques_supprimes']} supprimé(s), "
      f"{fasta_stats['doublons_mis_de_cote']} mis de côté")
print(f"  CSV   : {csv_stats_total['fichiers']} fichiers, "
      f"{csv_stats_total['lignes']} lignes uniques, "
      f"{csv_stats_total['doublons']} doublon(s)")
print(f"\n  Résultat → {OUT_DIR}/")
print(f"  Doublons → {DUP_DIR}/")
print("=" * 60)

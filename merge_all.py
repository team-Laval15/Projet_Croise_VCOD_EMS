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
for ds in CSV_DATASETS:
    ensure_dir(os.path.join(OUT_DIR, ds))
    ensure_dir(os.path.join(DUP_DIR, ds))

# ============================================================
# Résumé
# ============================================================

print("\n" + "=" * 60)
print("  Résumé")
print("=" * 60)
print(f"  FASTA anciens  : {FASTA_OLD_DIR}/")
print(f"  FASTA nouveaux : {FASTA_NEW_DIR}/")
print(f"  CSV anciens    : {CSV_OLD_DIR}/{{dataset}}/")
print(f"  CSV nouveaux   : {CSV_NEW_DIR}/{{dataset}}/")
print(f"  Sortie fusion  : {OUT_DIR}/")
print(f"  Doublons       : {DUP_DIR}/")
print("=" * 60)

confirm = input("\nLancer la fusion ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)

print()

# ============================================================
# Helpers
# ============================================================

def list_files(folder, extension):
    """Liste tous les fichiers avec l'extension donnée dans un dossier."""
    if not os.path.isdir(folder):
        return []
    return [f for f in os.listdir(folder) if f.endswith(extension)]


def read_fasta(path):
    """
    Lit un fichier FASTA et retourne une liste de (header, sequence).
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
                    sequences.append((header, seq))
                header = line
                seq    = ""
            else:
                seq += line

    if header is not None:
        sequences.append((header, seq))

    return sequences


def extract_fasta_id(header):
    """Extrait l'identifiant depuis un header FASTA (premier champ avant le point)."""
    # >U58801.23.7.1987.FR  →  U58801
    return header[1:].split(".")[0]


# ============================================================
# Fusion FASTA
# ============================================================

print("─" * 60)
print("  Fusion FASTA")
print("─" * 60)

old_fasta_files = set(list_files(FASTA_OLD_DIR, ".fasta"))
new_fasta_files = set(list_files(FASTA_NEW_DIR, ".fasta"))
all_fasta_files = old_fasta_files | new_fasta_files

fasta_stats = {"fichiers": 0, "sequences": 0, "doublons": 0}

for filename in sorted(all_fasta_files):
    old_path = os.path.join(FASTA_OLD_DIR, filename)
    new_path = os.path.join(FASTA_NEW_DIR, filename)

    sequences_vues = {}   # id → (header, seq)  — première occurrence gardée
    doublons       = {}   # id → (header, seq)   — occurrences en doublon

    # Lecture dans l'ordre : anciens d'abord, nouveaux ensuite
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
                doublons[seq_id] = (header, seq)

    # Écriture du fichier fusionné
    out_path = os.path.join(OUT_DIR, "sequences", filename)
    with open(out_path, "w", encoding="utf-8") as f:
        for seq_id, (header, seq) in sequences_vues.items():
            f.write(f"{header}\n{seq}\n")

    # Écriture des doublons
    if doublons:
        dup_path = os.path.join(DUP_DIR, "sequences", filename)
        with open(dup_path, "w", encoding="utf-8") as f:
            for seq_id, (header, seq) in doublons.items():
                f.write(f"{header}\n{seq}\n")

    fasta_stats["fichiers"]  += 1
    fasta_stats["sequences"] += len(sequences_vues)
    fasta_stats["doublons"]  += len(doublons)

    status = f"  ✔  {filename}  ({len(sequences_vues)} séquences"
    if doublons:
        status += f", {len(doublons)} doublon(s) mis de côté"
    print(status + ")")

print(f"\n  → {fasta_stats['fichiers']} fichier(s) fusionné(s), "
      f"{fasta_stats['sequences']} séquences uniques, "
      f"{fasta_stats['doublons']} doublon(s) au total")

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

    old_csv_files = set(list_files(old_ds_dir, ".csv"))
    new_csv_files = set(list_files(new_ds_dir, ".csv"))
    all_csv_files = old_csv_files | new_csv_files

    if not all_csv_files:
        print(f"  ⚠  Aucun fichier CSV trouvé pour [{dataset}]")
        continue

    ds_stats = {"fichiers": 0, "lignes": 0, "doublons": 0}

    print(f"\n  [{dataset}]")

    for filename in sorted(all_csv_files):
        old_path = os.path.join(old_ds_dir, filename)
        new_path = os.path.join(new_ds_dir, filename)

        rows_vues    = {}   # id → row dict  — première occurrence gardée
        doublons_csv = {}   # id → row dict  — occurrences en doublon
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
                        # Pas de colonne id → on garde tout sans dédup
                        fake_id = f"__noid_{len(rows_vues)}"
                        rows_vues[fake_id] = row
                    elif row_id not in rows_vues:
                        rows_vues[row_id] = row
                    else:
                        doublons_csv[row_id] = row

        if not fieldnames:
            print(f"    ⚠  {filename} — impossible de lire les colonnes, ignoré")
            continue

        # Écriture du fichier fusionné
        out_path = os.path.join(OUT_DIR, dataset, filename)
        with open(out_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows_vues.values())

        # Écriture des doublons
        if doublons_csv:
            dup_path = os.path.join(DUP_DIR, dataset, filename)
            with open(dup_path, "w", encoding="utf-8", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(doublons_csv.values())

        ds_stats["fichiers"]  += 1
        ds_stats["lignes"]    += len(rows_vues)
        ds_stats["doublons"]  += len(doublons_csv)

        status = f"    ✔  {filename}  ({len(rows_vues)} lignes"
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
      f"{fasta_stats['doublons']} doublon(s)")
print(f"  CSV   : {csv_stats_total['fichiers']} fichiers, "
      f"{csv_stats_total['lignes']} lignes uniques, "
      f"{csv_stats_total['doublons']} doublon(s)")
print(f"\n  Résultat → {OUT_DIR}/")
print(f"  Doublons → {DUP_DIR}/")
print("=" * 60)

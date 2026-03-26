# detect_amorces.py
# Détecte les amorces de la forme : G{4,10} + milieu(max 10 nt, max 3 G) + G{4,10}
#
# Règles précises sur le milieu :
#   - longueur 1 à 10 nucléotides
#   - commence par un nucléotide NON-G (sinon absorbé dans le bloc G initial)
#   - finit   par un nucléotide NON-G (sinon absorbé dans le bloc G final)
#   - contient au plus 3 'G' au total

import os
import re
import csv
from collections import defaultdict

# Le milieu doit commencer et finir par un non-G pour ne pas être
# absorbé dans les blocs G adjacents.
# Deux formes possibles :
#   - 1 seul caractère non-G          : [^g\n]
#   - 2+ caractères, bords non-G      : [^g\n][acgt]{0,8}[^g\n]
PATTERN = re.compile(
    r'(g{4,10})'
    r'((?:[^g\n][acgt]{0,8}[^g\n])|[^g\n])'
    r'(g{4,10})',
    re.IGNORECASE
)

COLONNES        = ["fichier", "id_sequence", "amorce", "position_debut", "position_fin"]
COLONNES_RECAP  = ["id_sequence", "nb_amorces"]


def lire_fasta(path):
    """Lit un fichier FASTA, retourne une liste de (id, sequence)."""
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
                    sequences.append((header[1:].split(".")[0], seq))
                header = line
                seq    = ""
            else:
                seq += line

    if header is not None:
        sequences.append((header[1:].split(".")[0], seq))

    return sequences


def detecter_amorces(seq_id, sequence):
    """
    Cherche toutes les amorces valides dans une séquence.
    Avance position par position pour capturer les chevauchements.
    """
    resultats = []
    seq_lower = sequence.lower()
    pos       = 0

    while pos < len(seq_lower):
        m = PATTERN.search(seq_lower, pos)
        if not m:
            break

        milieu = m.group(2)

        # Validation finale : au plus 3 G dans le milieu
        if milieu.lower().count("g") <= 3:
            resultats.append({
                "id_sequence":    seq_id,
                "amorce":         m.group(0),
                "position_debut": m.start() + 1,  # 1-indexé
                "position_fin":   m.end(),
            })

        pos = m.start() + 1

    return resultats


# ============================================================
# Paramètres
# ============================================================

print("=" * 60)
print("  Detection d'amorces dans les sequences FASTA")
print("=" * 60)
print("\nForme : G{4-10} + milieu(1-10 nt, debut/fin non-G, max 3 G) + G{4-10}")
print("Exemple : gggggacgtcacgcggggg\n")

FASTA_DIR   = input("Dossier FASTA  : ").strip()
OUTPUT_FILE = input("Fichier de sortie CSV (defaut: amorces.csv) : ").strip() or "amorces.csv"

recap_input = input("Exporter un recap du nombre d'amorces par individu ? (o/n) : ").strip().lower()
RECAP       = recap_input == "o"

if RECAP:
    default_recap = OUTPUT_FILE.replace(".csv", "_par_individu.csv")
    recap_input2  = input(f"Fichier recap par individu (defaut: {default_recap}) : ").strip()
    RECAP_FILE    = recap_input2 if recap_input2 else default_recap

if not os.path.isdir(FASTA_DIR):
    print(f"\nDossier introuvable : {FASTA_DIR}")
    exit(1)

fichiers = sorted(f for f in os.listdir(FASTA_DIR) if f.endswith(".fasta"))
if not fichiers:
    print(f"\nAucun fichier .fasta dans : {FASTA_DIR}")
    exit(1)

print(f"\n{len(fichiers)} fichier(s) FASTA trouves")
confirm = input("Lancer la detection ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annule.")
    exit(0)

# ============================================================
# Détection
# ============================================================

print()
stats             = {"fichiers": 0, "sequences": 0, "amorces": 0}
resultats_globaux = []
amorces_par_id    = defaultdict(int)   # { seq_id: nb_amorces }

for fname in fichiers:
    path      = os.path.join(FASTA_DIR, fname)
    sequences = lire_fasta(path)
    nb_fichier = 0

    for seq_id, sequence in sequences:
        stats["sequences"] += 1
        amorces = detecter_amorces(seq_id, sequence)

        for a in amorces:
            resultats_globaux.append({
                "fichier":        fname,
                "id_sequence":    a["id_sequence"],
                "amorce":         a["amorce"],
                "position_debut": a["position_debut"],
                "position_fin":   a["position_fin"],
            })
            amorces_par_id[seq_id] += 1
            nb_fichier += 1

    stats["fichiers"] += 1
    stats["amorces"]  += nb_fichier

    if nb_fichier > 0:
        print(f"  [OK] {fname} — {len(sequences)} sequences, {nb_fichier} amorce(s)")
    else:
        print(f"  [--] {fname} — {len(sequences)} sequences, aucune amorce")

# ============================================================
# Export CSV principal
# ============================================================

with open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=COLONNES)
    writer.writeheader()
    writer.writerows(resultats_globaux)

# ============================================================
# Export récap par individu (optionnel)
# ============================================================

if RECAP:
    # Trier par nb_amorces décroissant pour faciliter la lecture
    recap_rows = sorted(
        [{"id_sequence": sid, "nb_amorces": nb} for sid, nb in amorces_par_id.items()],
        key=lambda x: x["nb_amorces"],
        reverse=True
    )

    with open(RECAP_FILE, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLONNES_RECAP)
        writer.writeheader()
        writer.writerows(recap_rows)

    print(f"\n  Recap par individu : {RECAP_FILE}  ({len(recap_rows)} individu(s) avec amorces)")

# ============================================================
# Fin
# ============================================================

print(f"\n{'='*60}")
print(f"  Termine")
print(f"  Fichiers analyses    : {stats['fichiers']}")
print(f"  Sequences analysees  : {stats['sequences']}")
print(f"  Amorces trouvees     : {stats['amorces']}")
if RECAP:
    print(f"  Individus distincts  : {len(amorces_par_id)}")
print(f"  Resultats            : {OUTPUT_FILE}")
if RECAP:
    print(f"  Recap individus      : {RECAP_FILE}")
print(f"{'='*60}")

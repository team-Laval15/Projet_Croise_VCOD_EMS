# stream_pipeline.py

import os
import re
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pipelines.fetch_data  import fetch_fasta, fetch_csv, FASTA_DATASETS, CSV_DATASETS
from pipelines.fetch_data  import OK, EMPTY, TOO_EARLY, ERROR
from pipelines.split_fasta import split_fasta
from pipelines.utils       import ensure_dir

ALL_DATASETS        = FASTA_DATASETS | CSV_DATASETS
DEFAULT_RETRY_DELAY = 3 * 60
DEFAULT_WORKERS     = 10


# ============================================================
# Utilitaires
# ============================================================

def iter_trimesters(start_year, start_q, end_year, end_q):
    """Génère tous les (année, trimestre) dans l'ordre chronologique."""
    year, q = start_year, start_q
    while (year, q) <= (end_year, end_q):
        yield year, q
        q += 1
        if q > 4:
            q = 1
            year += 1


def detect_resume_point(base_dir, datasets, countries):
    """
    Parcourt les dossiers existants pour trouver le dernier trimestre
    déjà téléchargé (au moins un fichier présent pour ce trimestre).

    Structure attendue :
      - CSV  : base_dir/{dataset}/{country}/{dataset}_{country}_{year}Q{q}.csv
      - FASTA: base_dir/sequences/{country}_{year}Q{q}.fasta

    Retourne (last_year, last_q) ou None si rien trouvé.
    """
    pattern_csv   = re.compile(r'_(\d{4})Q([1-4])\.csv$')
    pattern_fasta = re.compile(r'_(\d{4})Q([1-4])\.fasta$')

    found = set()  # ensemble de (year, q)

    for dataset in datasets:
        ds_dir = os.path.join(base_dir, dataset)
        if not os.path.isdir(ds_dir):
            continue

        if dataset == "sequences":
            # Fichiers à plat dans sequences/
            for fname in os.listdir(ds_dir):
                m = pattern_fasta.search(fname)
                if m and os.path.getsize(os.path.join(ds_dir, fname)) > 0:
                    found.add((int(m.group(1)), int(m.group(2))))
        else:
            # Sous-dossiers par pays
            for country in countries:
                country_dir = os.path.join(ds_dir, country)
                if not os.path.isdir(country_dir):
                    continue
                for fname in os.listdir(country_dir):
                    m = pattern_csv.search(fname)
                    if m and os.path.getsize(os.path.join(country_dir, fname)) > 0:
                        found.add((int(m.group(1)), int(m.group(2))))

    if not found:
        return None

    return max(found)  # (year, q) le plus récent


def already_done(base_dir, datasets, country, year, quarter):
    """
    Retourne True si tous les datasets demandés ont déjà un fichier
    non-vide pour ce pays/trimestre.
    """
    trimester = f"{year}Q{quarter}"
    for dataset in datasets:
        if dataset == "sequences":
            path = os.path.join(base_dir, "sequences", f"{country}_{trimester}.fasta")
        else:
            path = os.path.join(base_dir, dataset, country,
                                f"{dataset}_{country}_{trimester}.csv")
        if not (os.path.isfile(path) and os.path.getsize(path) > 0):
            return False
    return True


# ============================================================
# Fetch
# ============================================================

def fetch_country_trimester(server, datasets, country, year, quarter,
                             dataset_dirs, token):
    """
    Fetch tous les datasets pour un pays + trimestre.
    Retourne (country, statuts_par_dataset, chemins_par_dataset).
    dataset_dirs : { dataset: dossier_racine_du_dataset }
    """
    statuts = {}
    chemins = {}

    for dataset in datasets:
        base_dir = dataset_dirs[dataset]
        if dataset == "sequences":
            status, path = fetch_fasta(
                server, country, year, quarter, base_dir, token=token
            )
        else:
            status, path = fetch_csv(
                server, dataset, country, year, quarter, base_dir, token=token
            )
        statuts[dataset] = status
        chemins[dataset] = path

    return country, statuts, chemins


def fetch_all_countries(server, datasets, countries, year, quarter,
                        dataset_dirs, token, max_workers):
    """
    Fetch un trimestre pour tous les pays en parallèle.
    Retourne all_results = { country: (statuts, chemins) }
    """
    all_results = {}

    with ThreadPoolExecutor(max_workers=min(max_workers, len(countries))) as executor:
        futures = {
            executor.submit(
                fetch_country_trimester,
                server, datasets, country, year, quarter,
                dataset_dirs, token
            ): country
            for country in countries
        }
        for future in as_completed(futures):
            country, statuts, chemins = future.result()
            all_results[country] = (statuts, chemins)

    return all_results


def trimester_ready(all_results, datasets):
    """
    Le trimestre est considéré prêt si aucun pays n'a renvoyé TOO_EARLY.
    """
    for country, (statuts, _) in all_results.items():
        for dataset in datasets:
            if statuts.get(dataset) == TOO_EARLY:
                return False
    return True


# ============================================================
# Paramètres
# ============================================================

print("=" * 60)
print("     Pipeline streaming — tous pays, temps réel")
print("=" * 60)

SERVER = input("\nServeur (ex: http://ton_serveur) : ").strip()
TOKEN  = input("Token d'accès (laisser vide si aucun) : ").strip() or None

print(f"\nDatasets disponibles : {', '.join(sorted(ALL_DATASETS))}")
datasets_input = input("Datasets à récupérer (séparés par des virgules) : ").strip().lower()
DATASETS = {d.strip() for d in datasets_input.split(",") if d.strip() in ALL_DATASETS}
invalid  = {d.strip() for d in datasets_input.split(",") if d.strip() not in ALL_DATASETS}
if invalid:
    print(f"⚠  Datasets inconnus ignorés : {', '.join(invalid)}")
if not DATASETS:
    print("❌ Aucun dataset valide. Arrêt.")
    exit(1)

print("\nEntrer les pays séparés par des virgules (ex: CM,CA,FR,DE,GB)")
COUNTRIES = [c.strip() for c in input("Pays : ").strip().upper().split(",") if c.strip()]
if not COUNTRIES:
    print("❌ Aucun pays renseigné. Arrêt.")
    exit(1)

START_YEAR = int(input("\nAnnée de début : ").strip())
START_Q    = int(input("Trimestre de début (1-4) : ").strip())
END_YEAR   = int(input("Année de fin : ").strip())
END_Q      = int(input("Trimestre de fin (1-4) : ").strip())

BASE_DIR = input("\nDossier racine de sortie (défaut: data_final) : ").strip() or "data_final"

retry_input   = input(f"\nDélai entre retries en minutes (défaut: 3)  : ").strip()
RETRY_DELAY   = int(retry_input) * 60 if retry_input else DEFAULT_RETRY_DELAY

workers_input = input(f"Requêtes parallèles max (défaut: 10)        : ").strip()
MAX_WORKERS   = int(workers_input) if workers_input else DEFAULT_WORKERS

# Un sous-dossier racine par dataset (les sous-dossiers pays sont créés dans fetch_csv)
DATASET_DIRS = {}
for ds in DATASETS:
    folder = os.path.join(BASE_DIR, ds)
    ensure_dir(folder)
    DATASET_DIRS[ds] = folder

# ============================================================
# Reprise automatique
# ============================================================

resume = detect_resume_point(BASE_DIR, DATASETS, COUNTRIES)
EFFECTIVE_START_YEAR = START_YEAR
EFFECTIVE_START_Q    = START_Q

if resume:
    last_year, last_q = resume
    # On repart du trimestre suivant le dernier trouvé
    next_q    = last_q + 1
    next_year = last_year
    if next_q > 4:
        next_q    = 1
        next_year += 1

    # On ne reprend que si le point de reprise est dans la plage demandée
    if (next_year, next_q) > (START_YEAR, START_Q):
        EFFECTIVE_START_YEAR = next_year
        EFFECTIVE_START_Q    = next_q
        print(f"\n♻  Reprise détectée : dernier trimestre trouvé = {last_year}Q{last_q}")
        print(f"   → Démarrage à partir de {EFFECTIVE_START_YEAR}Q{EFFECTIVE_START_Q}")
    else:
        print(f"\n♻  Reprise détectée ({last_year}Q{last_q}) mais déjà avant le début demandé — on repart du début.")

# ============================================================
# Résumé
# ============================================================

trimesters  = list(iter_trimesters(EFFECTIVE_START_YEAR, EFFECTIVE_START_Q, END_YEAR, END_Q))
nb_requetes = len(trimesters) * len(COUNTRIES) * len(DATASETS)
cout_estime = nb_requetes * 0.01

print("\n" + "=" * 60)
print("  Résumé")
print("=" * 60)
print(f"  Serveur      : {SERVER}")
print(f"  Token        : {'✔ défini' if TOKEN else '✘ absent'}")
print(f"  Datasets     : {', '.join(sorted(DATASETS))}")
print(f"  Pays         : {', '.join(COUNTRIES)} ({len(COUNTRIES)} pays)")
print(f"  Période      : {EFFECTIVE_START_YEAR}Q{EFFECTIVE_START_Q} → {END_YEAR}Q{END_Q}  ({len(trimesters)} trimestres)")
if (EFFECTIVE_START_YEAR, EFFECTIVE_START_Q) != (START_YEAR, START_Q):
    print(f"  (demandé dès  {START_YEAR}Q{START_Q}, reprise automatique appliquée)")
print(f"  Requêtes     : ~{nb_requetes}  (~{cout_estime:.2f}€ estimés, hors retries)")
print(f"  Parallélisme : {MAX_WORKERS} requêtes simultanées par trimestre")
print(f"  Retry        : toutes les {RETRY_DELAY // 60} min si 'too early'")
print(f"  Structure    :")
print(f"    CSV  → {BASE_DIR}/{{dataset}}/{{pays}}/{{dataset}}_{{pays}}_{{trimestre}}.csv")
print(f"    FASTA→ {BASE_DIR}/sequences/{{pays}}_{{trimestre}}.fasta  (à plat)")
print("=" * 60)

if not trimesters:
    print("\n✅ Tous les trimestres sont déjà téléchargés. Rien à faire.")
    exit(0)

confirm = input("\nLancer le pipeline ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)

# ============================================================
# Pipeline principal
# ============================================================

print(f"\n🚀 Démarrage — {datetime.now().strftime('%H:%M:%S')}\n")

done       = 0
all_errors = []   # [(label, country, dataset, message)]

for year, quarter in trimesters:
    label    = f"{year}Q{quarter}"
    attempts = 0

    while True:
        attempts += 1
        now = datetime.now().strftime("%H:%M:%S")
        print(f"[{now}] ⏳ {label} — tentative {attempts}  ({len(COUNTRIES)} pays × {len(DATASETS)} datasets)")

        all_results = fetch_all_countries(
            SERVER, DATASETS, COUNTRIES, year, quarter,
            DATASET_DIRS, TOKEN, MAX_WORKERS
        )

        if not trimester_ready(all_results, DATASETS):
            early_countries = [
                c for c, (statuts, _) in all_results.items()
                if any(s == TOO_EARLY for s in statuts.values())
            ]
            now = datetime.now().strftime("%H:%M:%S")
            print(f"[{now}] 💤 {label} pas encore dispo "
                  f"({len(early_countries)} pays en attente : {', '.join(sorted(early_countries))}) "
                  f"— retry dans {RETRY_DELAY // 60} min\n")
            time.sleep(RETRY_DELAY)
            continue

        # Trimestre prêt : on traite les résultats
        ok_countries    = []
        empty_countries = []
        error_countries = []

        for country, (statuts, chemins) in all_results.items():
            statut_list = list(statuts.values())

            # Pays entièrement vide → on ignore, on ne retente pas
            if all(s == EMPTY for s in statut_list):
                empty_countries.append(country)
                continue

            # Pays avec au moins une erreur → on log et on passe
            if any(s == ERROR for s in statut_list):
                error_countries.append(country)
                for ds, s in statuts.items():
                    if s == ERROR:
                        all_errors.append((label, country, ds, "ERROR retourné par fetch"))

            # Split FASTA si sequences récupéré avec succès
            if "sequences" in DATASETS and statuts.get("sequences") == OK and chemins.get("sequences"):
                split_fasta(chemins["sequences"], DATASET_DIRS["sequences"])

            if any(s == OK for s in statut_list):
                ok_countries.append(country)

        done += 1
        now = datetime.now().strftime("%H:%M:%S")
        print(f"[{now}] ✅ {label} traité ({done}/{len(trimesters)})")
        print(f"         ✔  {len(ok_countries)} pays avec données  : {', '.join(sorted(ok_countries))}")
        if empty_countries:
            print(f"         ○  {len(empty_countries)} pays sans données : {', '.join(sorted(empty_countries))}")
        if error_countries:
            print(f"         ✘  {len(error_countries)} pays en erreur   : {', '.join(sorted(error_countries))}")
        print()
        break

# ============================================================
# Fin
# ============================================================

print(f"{'='*60}")
print(f"  ✅ Pipeline terminé — {datetime.now().strftime('%H:%M:%S')}")
print(f"  Données disponibles dans : {BASE_DIR}/")
print(f"    CSV  → {BASE_DIR}/{{dataset}}/{{pays}}/")
print(f"    FASTA→ {BASE_DIR}/sequences/  (à plat)")
print(f"{'='*60}")

# Récap des erreurs
if all_errors:
    print(f"\n⚠  {len(all_errors)} erreur(s) rencontrée(s) (pays ignorés) :")
    print(f"  {'Trimestre':<10} {'Pays':<6} {'Dataset':<14} Raison")
    print(f"  {'-'*50}")
    for lbl, country, ds, msg in sorted(all_errors):
        print(f"  {lbl:<10} {country:<6} {ds:<14} {msg}")
    print()
else:
    print("\n  Aucune erreur. 🎉")

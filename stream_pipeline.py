# stream_pipeline.py

import os
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


def iter_trimesters(start_year, start_q, end_year, end_q):
    """Génère tous les (année, trimestre) dans l'ordre chronologique."""
    year, q = start_year, start_q
    while (year, q) <= (end_year, end_q):
        yield year, q
        q += 1
        if q > 4:
            q = 1
            year += 1


def fetch_country_trimester(server, datasets, country, year, quarter,
                             dataset_dirs, token):
    """
    Fetch tous les datasets pour un pays + trimestre.
    Retourne (country, statuts_par_dataset, chemins_par_dataset).
    dataset_dirs : { dataset: dossier_de_sortie }
    """
    statuts = {}
    chemins = {}

    for dataset in datasets:
        output_dir = dataset_dirs[dataset]
        if dataset == "sequences":
            status, path = fetch_fasta(
                server, country, year, quarter, output_dir, token=token
            )
        else:
            status, path = fetch_csv(
                server, dataset, country, year, quarter, output_dir, token=token
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
    Les pays EMPTY ou ERROR sont acceptés (pas de données = données absentes, on continue).
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

# Un sous-dossier par dataset
DATASET_DIRS = {}
for ds in DATASETS:
    folder = os.path.join(BASE_DIR, ds)
    ensure_dir(folder)
    DATASET_DIRS[ds] = folder

# ============================================================
# Résumé
# ============================================================

trimesters  = list(iter_trimesters(START_YEAR, START_Q, END_YEAR, END_Q))
nb_requetes = len(trimesters) * len(COUNTRIES) * len(DATASETS)
cout_estime = nb_requetes * 0.01

print("\n" + "=" * 60)
print("  Résumé")
print("=" * 60)
print(f"  Serveur      : {SERVER}")
print(f"  Token        : {'✔ défini' if TOKEN else '✘ absent'}")
print(f"  Datasets     : {', '.join(sorted(DATASETS))}")
print(f"  Pays         : {', '.join(COUNTRIES)} ({len(COUNTRIES)} pays)")
print(f"  Période      : {START_YEAR}Q{START_Q} → {END_YEAR}Q{END_Q}  ({len(trimesters)} trimestres)")
print(f"  Requêtes     : ~{nb_requetes}  (~{cout_estime:.2f}€ estimés, hors retries)")
print(f"  Parallélisme : {MAX_WORKERS} requêtes simultanées par trimestre")
print(f"  Retry        : toutes les {RETRY_DELAY // 60} min si 'too early'")
print(f"  Dossiers     :")
for ds, folder in sorted(DATASET_DIRS.items()):
    print(f"    {ds:12} → {folder}/")
print("=" * 60)

confirm = input("\nLancer le pipeline ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)

# ============================================================
# Pipeline principal
# ============================================================

print(f"\n🚀 Démarrage — {datetime.now().strftime('%H:%M:%S')}\n")

done = 0

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
            # Au moins un pays a répondu "too early" → on attend et on réessaie
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
            # Détermine le statut global du pays sur ce trimestre
            statut_list = list(statuts.values())
            if all(s == EMPTY for s in statut_list):
                empty_countries.append(country)
                continue
            if any(s == ERROR for s in statut_list):
                error_countries.append(country)

            # Split FASTA si sequences récupéré
            if "sequences" in DATASETS and chemins.get("sequences"):
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
for ds, folder in sorted(DATASET_DIRS.items()):
    print(f"    {ds:12} → {folder}/")
print(f"{'='*60}")

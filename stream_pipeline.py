# stream_pipeline.py
 
import time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pipelines.fetch_data  import fetch_fasta, fetch_csv, FASTA_DATASETS, CSV_DATASETS
from pipelines.split_fasta import split_fasta
from pipelines.merge_fasta import merge_directories
from pipelines.utils       import ensure_dir
 
ALL_DATASETS        = FASTA_DATASETS | CSV_DATASETS
DEFAULT_RETRY_DELAY = 3 * 60  # secondes
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
                             temp_dir, csv_dir, token):
    """
    Fetch tous les datasets pour un pays + trimestre.
    Retourne (country, success, results_dict).
    """
    results = {}
    for dataset in datasets:
        if dataset == "sequences":
            path = fetch_fasta(server, country, year, quarter, temp_dir, token=token)
        else:
            path = fetch_csv(server, dataset, country, year, quarter, csv_dir, token=token)
        results[dataset] = path
 
    success = all(p is not None for p in results.values())
    return country, success, results
 
 
def fetch_all_countries(server, datasets, countries, year, quarter,
                        temp_dir, csv_dir, token, max_workers):
    """
    Fetch un trimestre pour tous les pays en parallèle.
    Retourne (tous_ok, results_par_pays).
    """
    all_results = {}
 
    with ThreadPoolExecutor(max_workers=min(max_workers, len(countries))) as executor:
        futures = {
            executor.submit(
                fetch_country_trimester,
                server, datasets, country, year, quarter,
                temp_dir, csv_dir, token
            ): country
            for country in countries
        }
        for future in as_completed(futures):
            country, success, results = future.result()
            all_results[country] = (success, results)
 
    tous_ok = all(s for s, _ in all_results.values())
    return tous_ok, all_results
 
 
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
 
TEMP_DIR     = input("\nDossier temporaire FASTA bruts  : ").strip() or "new_data"
SPLIT_DIR    = input("Dossier FASTA splités           : ").strip() or "fasta_split"
MERGE_DIR    = input("Dossier FASTA final (merge)     : ").strip() or "fasta_merge"
OLD_DATA_DIR = input("Dossier anciennes données FASTA : ").strip() or "fasta_data"
CSV_DIR      = input("Dossier CSV                     : ").strip() or "csv_data"
 
retry_input   = input(f"\nDélai entre retries en minutes (défaut: 3)  : ").strip()
RETRY_DELAY   = int(retry_input) * 60 if retry_input else DEFAULT_RETRY_DELAY
 
workers_input = input(f"Requêtes parallèles max (défaut: 10)        : ").strip()
MAX_WORKERS   = int(workers_input) if workers_input else DEFAULT_WORKERS
 
for d in (TEMP_DIR, SPLIT_DIR, MERGE_DIR, CSV_DIR):
    ensure_dir(d)
 
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
print(f"  Retry        : toutes les {RETRY_DELAY // 60} min si trimestre pas encore dispo")
print("=" * 60)
 
confirm = input("\nLancer le pipeline ? (o/n) : ").strip().lower()
if confirm != "o":
    print("Annulé.")
    exit(0)
 
# ============================================================
# Pipeline principal
# ============================================================
 
print(f"\n🚀 Démarrage — {datetime.now().strftime('%H:%M:%S')}\n")
 
done          = 0
fasta_fetched = []
 
for year, quarter in trimesters:
    label    = f"{year}Q{quarter}"
    attempts = 0
 
    while True:
        attempts += 1
        now = datetime.now().strftime("%H:%M:%S")
        print(f"[{now}] ⏳ {label} — tentative {attempts}  ({len(COUNTRIES)} pays en parallèle)")
 
        tous_ok, all_results = fetch_all_countries(
            SERVER, DATASETS, COUNTRIES, year, quarter,
            TEMP_DIR, CSV_DIR, TOKEN, MAX_WORKERS
        )
 
        if tous_ok:
            # Split FASTA pour chaque pays
            if "sequences" in DATASETS:
                for country, (_, results) in all_results.items():
                    fasta_path = results.get("sequences")
                    if fasta_path:
                        split_fasta(fasta_path, SPLIT_DIR)
                        fasta_fetched.append(fasta_path)
 
            done += 1
            now            = datetime.now().strftime("%H:%M:%S")
            ok_list        = sorted(c for c, (s, _) in all_results.items() if s)
            fail_list      = sorted(c for c, (s, _) in all_results.items() if not s)
 
            print(f"[{now}] ✅ {label} terminé ({done}/{len(trimesters)})")
            print(f"         ✔  {len(ok_list)} pays OK : {', '.join(ok_list)}")
            if fail_list:
                print(f"         ✘  {len(fail_list)} pays en échec : {', '.join(fail_list)}")
            print()
            break
 
        else:
            fail_list = sorted(c for c, (s, _) in all_results.items() if not s)
            now       = datetime.now().strftime("%H:%M:%S")
            print(f"[{now}] 💤 {label} pas encore dispo "
                  f"({len(fail_list)} pays manquants : {', '.join(fail_list)}) "
                  f"— retry dans {RETRY_DELAY // 60} min\n")
            time.sleep(RETRY_DELAY)
 
# ============================================================
# Merge final FASTA
# ============================================================
 
if "sequences" in DATASETS and fasta_fetched:
    print(f"\n{'='*60}")
    print("  Fusion finale des fichiers FASTA")
    print(f"{'='*60}")
    merge_directories(OLD_DATA_DIR, SPLIT_DIR, MERGE_DIR)
 
print(f"\n{'='*60}")
print(f"  ✅ Pipeline terminé — {datetime.now().strftime('%H:%M:%S')}")
if "sequences" in DATASETS:
    print(f"  FASTA fusionnés → {MERGE_DIR}")
if DATASETS - {"sequences"}:
    print(f"  CSV             → {CSV_DIR}")
print(f"{'='*60}")

# Genome annotation pipeline

Nextflow workflow for Bakta annotation, geNomad virus/plasmid detection, MOB-suite plasmid reconstruction, StarAMR, RGI and AMRFinderPlus. Tools run in version-pinned containers; databases are prepared separately.

**Status: The pipeline completed a four-genome SLURM run on Nextflow 26.10.6 and SingularityCE 3.11.1 using Nextflow's legacy parser. Strict-parser compatibility after the recent syntax fix has not yet been checked. Local execution, empty-input handling and broader output validation also remain to be tested. Review results before interpreting them.

## Requirements

- Linux, Java 17+, Nextflow, SingularityCE/Apptainer; SLURM for cluster runs.
- Internet access for image and database downloads.
- Writable storage for images, databases and Nextflow work files. Database and work paths must be accessible from compute nodes.

## Install

Clone the repository and enter it:

```bash
git clone https://github.com/kayanac/bacterial-genome-annotation.git
cd bacterial-genome-annotation
export PROJECT_DIR="$PWD"
mkdir -p images
```

The workflow uses pinned registry images by default, so supplying local container paths is optional. Nextflow/Singularity can fetch the pinned images when needed. To download them in advance (for example, when compute nodes cannot access the registry), set writable cache and temporary paths and run:

```bash
export SINGULARITY_CACHEDIR=/path/to/singularity-cache
export SINGULARITY_TMPDIR=/path/to/singularity-build-tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

tail -n +2 containers.tsv | while IFS=$'\t' read -r tool version image digest; do
  singularity pull "images/${tool}.sif" "docker://${image%:*}@${digest}" </dev/null || exit 1
done
```

## Prepare databases

Choose a new database directory on writable storage accessible from compute nodes. The full Bakta database download is large (about 32 GB compressed in the tested setup), and extraction needs additional space. These commands download available database releases, which can change over time; keep the resulting database directories for reruns.

```bash
export DB_ROOT=/path/to/genome-annotation-databases
mkdir "$DB_ROOT" || { echo "Choose a new, unused DB_ROOT" >&2; }
export DB_ROOT="$(cd "$DB_ROOT" && pwd -P)"

export BAKTA_DB="$DB_ROOT/bakta/db"
export GENOMAD_DB="$DB_ROOT/genomad_db"
export MOB_DB="$DB_ROOT/mob-suite"
export STARAMR_DB="$DB_ROOT/staramr"
export AMRFINDER_ROOT="$DB_ROOT/amrfinder"
mkdir "$DB_ROOT/bakta" "$MOB_DB" "$AMRFINDER_ROOT"
```

Containerized database setup needs writable temporary storage visible inside the containers. Choose a suitable path for your machine or cluster; network filesystems may not support all temporary build operations.

```bash
export SETUP_TMP=/path/to/writable/setup-tmp
mkdir -p "$SETUP_TMP"
export SETUP_TMP="$(cd "$SETUP_TMP" && pwd -P)"
```

Run each setup command separately. Stop and resolve any error before continuing:

```bash
singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/Bakta.sif" \
  bakta_db download --output "$DB_ROOT/bakta" --type full

singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/geNomad.sif" genomad download-database "$DB_ROOT"

mkdir "$DB_ROOT/mob-setup-work"
(cd "$DB_ROOT/mob-setup-work" && singularity exec --cleanenv \
  --bind "$DB_ROOT:$DB_ROOT" --bind "$SETUP_TMP:$SETUP_TMP" \
  --env "TMPDIR=$SETUP_TMP" "$PROJECT_DIR/images/MOB-suite.sif" \
  mob_init --database_directory "$MOB_DB")

singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/StarAMR.sif" staramr db build --dir "$STARAMR_DB"

singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/AMRFinderPlus.sif" \
  amrfinder_update --database "$AMRFINDER_ROOT"
```

Resolve the standalone AMRFinderPlus database to its dated directory. Check database versions and keep the database snapshots unchanged for reproducible reruns. RGI's image includes its CARD database.

```bash
export AMRFINDER_DB="$(readlink -f "$AMRFINDER_ROOT/latest")"
test -f "$AMRFINDER_DB/version.txt" || echo "AMRFinderPlus database missing"

singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/StarAMR.sif" staramr db info "$STARAMR_DB"
singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  "$PROJECT_DIR/images/AMRFinderPlus.sif" \
  amrfinder --database "$AMRFINDER_DB" --database_version
singularity exec --cleanenv "$PROJECT_DIR/images/RGI.sif" rgi database --version
```

In the tested setup, Bakta used database 6.0 full, geNomad 1.9, AMRFinderPlus 2026-10-29.1 and RGI's bundled CARD 4.0.1. StarAMR may report non-default database revisions; review and record the `db info` output. Do not update a database snapshot used by an earlier run.

## Run on SLURM

Put one uncompressed `.fa`, `.fna` or `.fasta` assembly per sample in an input directory. Sample basenames must be unique and start with a letter or digit.

Run Nextflow from the project directory on the login node, following your cluster's policy. SLURM executes the workflow tasks. The database directories themselves are required for this full workflow, but their command-line path options are optional: the config reads these environment variables, or you can pass the paths as `--bakta_db`, `--genomad_db`, `--mob_db`, `--staramr_db` and `--amrfinder_db` on the run command.

```bash
export BAKTA_DB=/path/to/genome-annotation-databases/bakta/db
export GENOMAD_DB=/path/to/genome-annotation-databases/genomad_db
export MOB_DB=/path/to/genome-annotation-databases/mob-suite
export STARAMR_DB=/path/to/genome-annotation-databases/staramr
export AMRFINDER_DB=/path/to/genome-annotation-databases/amrfinder/DATE.VERSION

mkdir -p results
NXF_SYNTAX_PARSER=v1 nextflow -c nextflow_annotation.config run genome_annotation_ID_Cov.nf \
  -profile slurm \
  --input_assemblies /path/to/assemblies \
  --outDir "$PROJECT_DIR/results" \
  --run_name example_run \
  --threads 4 --cpus 4 --memory '32 GB' --time '6h' \
  -work-dir "$PROJECT_DIR/work-example"
```

The container path options (`--bakta_container`, `--genomad_container`, `--mob_container`, `--staramr_container`, `--rgi_container`, and `--amrfinder_container`) are also optional. The workflow defaults to the pinned registry images. Use those options only to override the defaults, such as when running a local `.sif` file; downloaded `.sif` files are not selected automatically.

`--samples sample1,sample2` optionally selects basenames without extensions. `--min_assembly_bases` defaults to `0`; a positive value skips shorter non-empty assemblies.

Keep the work directory so a failed run can be resumed with the same command plus `-resume`. Confirm Nextflow marks completed tasks as cached. Scheduler, container startup or staging errors can stop a run before the final status summary is produced.

## Results

Results are written under `results/<run_name>/`, including per-sample tool outputs, an ARG summary and `<run_name>_tool_status_summary.tsv`. Check the status summary:

```bash
cat "$PROJECT_DIR/results/example_run/example_run_tool_status_summary.tsv"
```

Tool failures and invalid empty assemblies should make the run exit nonzero. A no-hit result is not proof of biological absence. `MOB_NO_TYPER_RESULT` means MOB-suite produced no typing report; inspect its contig report and logs rather than interpreting it as proof that plasmids are absent.

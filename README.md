# Genome annotation pipeline

Nextflow workflow for Bakta annotation, geNomad virus/plasmid detection, MOB-suite plasmid reconstruction, StarAMR, RGI and AMRFinderPlus. Software runs in pinned containers; reference databases are stored separately.

**Status: draft.** A four-genome SLURM run completed with SingularityCE 3.11.1 and Nextflow 26.10.6 using the legacy parser. Results and scientific accuracy have not been comprehensively validated; local/Docker profiles and empty-input handling also need testing. Inspect outputs before using them for research.

## 1. Requirements

- Linux, Java 17+, Nextflow and SingularityCE; SLURM for the cluster profile. Docker is an alternative local profile, but its database setup is not covered below.
- Internet access for image and database downloads, and substantial free space. The tested full Bakta archive alone was about 32 GB compressed; extraction requires additional space.
- Shared project, database and work directories accessible on the login and compute nodes. Follow your cluster's policies on downloading and resource-intensive setup.

No separate Conda installation of the analysis tools is required. `containers.tsv` lists the pinned image tags and digests. `baseline_versions.tsv` is a historical migration reference, **not** the versions used by this workflow.

## 2. Get the project and images

Clone this repository using your GitHub authentication:

```bash
git clone https://github.com/kayanac/genome-annotation.git
cd genome-annotation
export PROJECT_DIR="$PWD"
mkdir -p images
```

Singularity needs space for its cache and temporary image builds. Choose writable locations suitable for your system. Temporary build storage must support the filesystem operations Singularity requires; network storage may not. For SLURM, use site-approved scratch and allocations where required.

```bash
# Replace these with suitable absolute paths.
export SINGULARITY_CACHEDIR=/path/to/singularity-cache
export SINGULARITY_TMPDIR=/path/to/singularity-build-tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"

tail -n +2 containers.tsv | while IFS=$'\t' read -r tool version image digest; do
  singularity pull "images/${tool}.sif" "docker://${image%:*}@${digest}" </dev/null || exit 1
done
ls -lh images/*.sif
```

The tag is removed from the pull reference because SingularityCE 3.11 rejects references containing both a tag and a digest. The digest remains pinned. Do not continue if a pull fails; resolve the error and pull the missing image.

## 3. Prepare databases

Use a **new** database directory and preserve it after setup. These commands download available releases, which may change over time. Do not overwrite an existing snapshot. Set up in a suitable SLURM allocation if required; `mob_init` can use up to 32 CPUs, limited by available CPU affinity.

```bash
# Replace with an unused path on a filesystem visible to all compute nodes.
export DB_ROOT=/path/to/genome-annotation-databases
mkdir "$DB_ROOT" || { echo "Choose a new, writable DB_ROOT" >&2; }
```

If `mkdir` fails, stop and choose another path. Then:

```bash
export DB_ROOT="$(cd "$DB_ROOT" && pwd -P)"
export BAKTA_DB="$DB_ROOT/bakta/db"
export GENOMAD_DB="$DB_ROOT/genomad_db"
export MOB_DB="$DB_ROOT/mob-suite"
export STARAMR_DB="$DB_ROOT/staramr"
export AMRFINDER_ROOT="$DB_ROOT/amrfinder"
mkdir "$DB_ROOT/bakta" "$MOB_DB" "$AMRFINDER_ROOT"

# Choose a writable temporary directory accessible from inside the container.
export SETUP_TMP=/path/to/writable/setup-tmp
mkdir -p "$SETUP_TMP"
export SETUP_TMP="$(cd "$SETUP_TMP" && pwd -P)"
```

Run each command separately and stop on error. Explicitly binding temporary storage is important when it is outside Singularity's automatic mounts.

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

StarAMR's destination must not already exist. Its build can warn that downloaded databases are non-default: check compatibility and record the revisions. If Bakta's internal AMRFinderPlus update fails after the main database download, inspect the error and retry the update rather than downloading Bakta again:

```bash
singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/Bakta.sif" \
  amrfinder_update --database "$BAKTA_DB/amrfinderplus-db"
```

Record versions and resolve the standalone AMRFinderPlus `latest` symlink to its dated directory:

```bash
export AMRFINDER_DB="$(readlink -f "$AMRFINDER_ROOT/latest")"
test -f "$AMRFINDER_DB/version.txt" || echo "AMRFinderPlus database missing"
cat "$BAKTA_DB/version.json" "$MOB_DB/status.txt"
singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  --bind "$SETUP_TMP:$SETUP_TMP" --env "TMPDIR=$SETUP_TMP" \
  "$PROJECT_DIR/images/StarAMR.sif" staramr db info "$STARAMR_DB"
singularity exec --cleanenv --bind "$DB_ROOT:$DB_ROOT" \
  "$PROJECT_DIR/images/AMRFinderPlus.sif" amrfinder --database "$AMRFINDER_DB" --database_version
singularity exec --cleanenv "$PROJECT_DIR/images/RGI.sif" rgi database --version
```

The tested setup used Bakta 6.0 full, geNomad 1.9, AMRFinderPlus 2026-10-29.1 and RGI's bundled CARD 4.0.1. StarAMR used non-default database revisions; check its `db info` and preserve the exact snapshots. A fresh setup may obtain different database contents, including Bakta's internal AMRFinderPlus database and MOB-suite taxonomy. For reproducible reruns preserve the entire database directories and symlinks, record versions and source revisions, and checksum files:

```bash
(cd "$DB_ROOT" && find bakta genomad_db mob-suite staramr amrfinder -type f -print0 \
  | sort -z | xargs -0 sha256sum) > "${DB_ROOT}.sha256sums.txt"
(cd "$DB_ROOT" && sha256sum -c "${DB_ROOT}.sha256sums.txt")
```

The checksum list covers regular files, not symlink targets; preserve and record symlinks separately. Database files and genomes should not be committed to GitHub.

## 4. Run

Place uncompressed `.fa`, `.fna`, or `.fasta` assemblies in an input directory. One file represents one sample; basenames must be unique and contain only letters, digits, dots, underscores or hyphens, starting with a letter or digit.

From the project directory, with the five database environment variables above exported:

```bash
mkdir -p results
NXF_SYNTAX_PARSER=v1 nextflow -c nextflow_annotation.config run genome_annotation_ID_Cov.nf \
  -profile slurm \
  --input_assemblies /path/to/assemblies \
  --outDir "$PROJECT_DIR/results" \
  --run_name example_run \
  --threads 4 --cpus 4 --memory '32 GB' --time '6h' \
  --bakta_container "$PROJECT_DIR/images/Bakta.sif" \
  --genomad_container "$PROJECT_DIR/images/geNomad.sif" \
  --mob_container "$PROJECT_DIR/images/MOB-suite.sif" \
  --staramr_container "$PROJECT_DIR/images/StarAMR.sif" \
  --rgi_container "$PROJECT_DIR/images/RGI.sif" \
  --amrfinder_container "$PROJECT_DIR/images/AMRFinderPlus.sif" \
  -work-dir "$PROJECT_DIR/work-example"
```

Use the login node to launch Nextflow under site policy; SLURM runs the tasks. The `v1` setting is needed for the tested Nextflow 26.10.6 because the workflow is not yet compatible with the default strict parser. Nextflow may warn about repeated parameter assignments. The local Singularity profile is `-profile local`; `docker_local` is available but untested. All database paths must be readable absolute paths visible on compute nodes and are mounted read-only.

By default all matching files are used. Select samples with `--samples sample1,sample2`, using basenames without extensions. `--min_assembly_bases` defaults to `0`; a positive value skips shorter non-empty assemblies. Keep the work directory and Nextflow cache; after resolving a failure, rerun the same command with `-resume`. Confirm tasks are marked cached before assuming they were reused.

## 5. Check results

Results are published under `results/<run_name>/`, including per-sample output directories, an ARG summary and `<run_name>_tool_status_summary.tsv`:

```bash
cat "$PROJECT_DIR/results/example_run/example_run_tool_status_summary.tsv"
```

The summary identifies the sample and tool. A tool command failure or invalid empty assembly should produce a nonzero final status; the other analyses are allowed to finish. Scheduler, container startup or staging errors can stop the run before the summary is written. Check Nextflow logs and per-task logs in the work directory.

`MOB_NO_TYPER_RESULT` means no typing report was generated, not proof that plasmids are absent; inspect the contig report and logs. A no-hit result is not proof of biological absence. Validate the ARG summary and individual tool outputs before interpreting them.

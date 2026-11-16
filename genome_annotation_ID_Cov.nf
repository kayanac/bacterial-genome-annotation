#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
  genome_annotation.nf

*/

params.input_assemblies = ''
params.outDir           = ''
params.run_name         = ''
params.samples          = ''     // OPTIONAL: sample name or comma-separated list (basename without extension)
params.threads          = 16
params.memory           = '64 GB'
params.time             = '12h'
params.min_assembly_bases = 0

params.cpus = params.cpus ?: params.threads
params.bakta_container = 'quay.io/biocontainers/bakta@sha256:d92c9fcf5e502b6607fb3672df98094ce98a0352ce98e72017f389db863d48a3'
params.genomad_container = 'quay.io/biocontainers/genomad@sha256:c648807baea891177fc571eab444a07ab310ce56c6f80b2939cc29ee892fa795'
params.mob_container = 'quay.io/biocontainers/mob_suite@sha256:0ae9dc0ee64c21e00e184e4bc02e20bc785a69382fed44ab69eb6c37aff5a69b'
params.staramr_container = 'quay.io/biocontainers/staramr@sha256:d3fd72ca3d635fd3563061f35760b7a6a016e6b801b14a99db1067e10c7da272'
params.rgi_container = 'quay.io/biocontainers/rgi@sha256:71d87a905d29bbc2dd2a44500bf99d5da2df4ccbef08123b047ad48350868583'
params.amrfinder_container = 'quay.io/biocontainers/ncbi-amrfinderplus@sha256:ea12f3f775905e6653ba089c46d471152e458c62733c0ca3677ced3ccf84e2be'

params.bakta_db = params.bakta_db ?: ''
params.genomad_db = params.genomad_db ?: ''
params.mob_db = params.mob_db ?: ''
params.staramr_db = params.staramr_db ?: ''
params.amrfinder_db = params.amrfinder_db ?: ''
/* Basic parameter validation */
if (!params.input_assemblies) {
    log.error "Missing required parameter: --input_assemblies (folder containing fasta files)"
    exit 1
}
if (!params.outDir) {
    log.error "Missing required parameter: --outDir (top output folder e.g. KADIR_GENOMIC_annotation_bysample)"
    exit 1
}
if (!params.run_name) {
    log.error "Missing required parameter: --run_name (subfolder under top output folder, e.g. 20250529)"
    exit 1
}
if (!(params.run_name ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
    log.error "--run_name may contain only letters, digits, dot, underscore, and hyphen"
    exit 1
}
try {
    params.threads = params.threads as Integer
    params.cpus = params.cpus as Integer
} catch (Exception e) {
    log.error "--threads and --cpus must be positive integers"
    exit 1
}
if (params.threads < 1 || params.cpus < 1) {
    log.error "--threads and --cpus must be positive integers"
    exit 1
}
try {
    params.min_assembly_bases = params.min_assembly_bases as Integer
} catch (Exception e) {
    log.error "--min_assembly_bases must be a non-negative integer"
    exit 1
}
if (params.min_assembly_bases < 0) {
    log.error "--min_assembly_bases must be a non-negative integer"
    exit 1
}
[
    bakta_db: params.bakta_db,
    genomad_db: params.genomad_db,
    mob_db: params.mob_db,
    staramr_db: params.staramr_db,
    amrfinder_db: params.amrfinder_db
].each { name, value ->
    if (!value) {
        log.error "Missing required database directory parameter: --${name}"
        exit 1
    }
    def database = file(value)
    if (!database.exists() || !database.isDirectory() || !database.toFile().canRead()) {
        log.error "Database directory must exist and be readable: --${name}=${value}"
        exit 1
    }
    if (!database.isAbsolute()) {
        log.error "Database paths must be absolute so container runtimes can bind them: --${name}=${value}"
        exit 1
    }
}
def asmDir = file(params.input_assemblies)
if (!asmDir.exists() || !asmDir.isDirectory() || !asmDir.toFile().canRead()) {
    log.error "Input assemblies folder must exist and be readable: ${params.input_assemblies}"
    exit 1
}

/* Create run root folder */
def run_root = file("${params.outDir}/${params.run_name}")
if (!run_root.exists() && !run_root.mkdirs()) {
    log.error "Could not create output directory: ${run_root}"
    exit 1
}
if (!run_root.isDirectory() || !run_root.toFile().canWrite()) {
    log.error "Output path must be a writable directory: ${run_root}"
    exit 1
}

/* Build channel of assemblies (case-insensitive extensions) */
def raw_assemblies_ch = Channel
    .fromPath("${params.input_assemblies}/*", type: 'file', checkIfExists: false)
    .filter { path -> path && path.toFile().canRead() && path.getName() ==~ /(?i).+\.(fa|fna|fasta)$/ }
    .map { path ->
        def nm = path.getName().replaceAll(/(?i)\.(fa|fna|fasta)$/, '')
        tuple(nm, path)
    }

/* If user supplied --samples, filter the channel to only those names */
if (params.samples?.trim()) {
    def wanted = params.samples.split(',').collect { it.trim() }
    raw_assemblies_ch = raw_assemblies_ch.filter { t -> wanted.contains(t[0]) }
}

/* materialize final channel and validate there are assemblies */
raw_assemblies_ch = raw_assemblies_ch
    .toList()
    .map { list ->
        if (!list || list.size() == 0) {
            throw new IllegalArgumentException("No readable FASTA files found in ${params.input_assemblies} after applying sample filter (${params.samples})")
        }
        def invalid = list.find { item -> !(item[0] ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/) }
        if (invalid) {
            throw new IllegalArgumentException("Sample names must start with a letter or digit and contain only letters, digits, dot, underscore, or hyphen: ${invalid[0]}")
        }
        def duplicate = list.groupBy { it[0] }.find { sample, items -> items.size() > 1 }
        if (duplicate) {
            throw new IllegalArgumentException("Multiple assembly files resolve to the same sample name: ${duplicate.key}")
        }
        return list
    }
    .flatMap { it }

raw_assemblies_ch.set { assemblies_ch }

/* ----------------------
   Bakta per-sample process
   ---------------------- */
process BAKTA_ANNOTATE {
    tag { sample_name }
    container params.bakta_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      tuple val(sample_name), path("${sample_name}/${sample_name}_bakta_annotation"), emit: bakta_out
      path "${sample_name}/${sample_name}_bakta_annotation/bakta_status.txt", emit: bakta_status

    script:
    """
    set -euo pipefail
    echo "BAKTA: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_bakta_annotation"
    TMP_PARENT=\$(mktemp -d "${sample_name}_bakta_parent.XXXXXX")
    TEMP_OUT="\$TMP_PARENT/${sample_name}_bakta_out"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/bakta_tmp.info"
    echo "TEMP_OUT=\$TEMP_OUT" >> "\$TMP_PARENT/bakta_tmp.info"

    echo "assembly: ${assembly}" > "\$TMP_PARENT/metadata.txt"
    echo "sample: ${sample_name}" >> "\$TMP_PARENT/metadata.txt"
    echo "bakta_db: ${params.bakta_db}" >> "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.bakta_container}" >> "\$TMP_PARENT/metadata.txt"
    echo "threads: ${task.cpus}" >> "\$TMP_PARENT/metadata.txt"

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "${assembly}" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    echo "seq_bases=\$seq_bases" >> "\$TMP_PARENT/metadata.txt"

    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/bakta_status.txt"
      echo ">${sample_name}_placeholder_annotation" > "\$TARGET_DIR/${sample_name}_placeholder.fasta"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/bakta_status.txt"
      echo ">${sample_name}_placeholder_annotation" > "\$TARGET_DIR/${sample_name}_placeholder.fasta"
      exit 0
    fi

    bakta --db "${params.bakta_db}" --output "\$TEMP_OUT" --threads ${task.cpus} --prefix "${sample_name}" "${assembly}" > "\$TMP_PARENT/bakta.stdout" 2> "\$TMP_PARENT/bakta.stderr" || rc=\$?

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      mkdir -p "\$TARGET_DIR"
      echo "BAKTA_FAILED rc=\$rc" > "\$TARGET_DIR/bakta_status.txt"
      echo "See tmp stderr at: \$TMP_PARENT/bakta.stderr" >> "\$TARGET_DIR/bakta_status.txt"
      cp -a "\$TMP_PARENT/bakta.stdout" "\$TMP_PARENT/bakta.stderr" "\$TARGET_DIR/" || true
      if [ -d "\$TEMP_OUT" ]; then
        cp -a "\$TEMP_OUT/." "\$TARGET_DIR/" || true
      fi
      exit 0
    fi

    if [ -d "\$TEMP_OUT" ]; then
      rm -rf "\$TARGET_DIR" || true
      mv "\$TEMP_OUT" "\$TARGET_DIR"
    else
      mkdir -p "\$TARGET_DIR"
      echo "BAKTA_FAILED: output directory not created: \$TEMP_OUT" > "\$TARGET_DIR/bakta_status.txt"
    fi

    if [ ! -s "\$TARGET_DIR/bakta_status.txt" ]; then
      echo "BAKTA_OK" > "\$TARGET_DIR/bakta_status.txt"
    fi
    """
}

/* ----------------------
   GENOMAD per-sample process 
   ---------------------- */
process GENOMAD_RUN {
    tag { sample_name }
    container params.genomad_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      path "${sample_name}/${sample_name}_genomad/${sample_name}_plasmid.fna", optional: true, emit: genomad_plasmid_fna
      path "${sample_name}/${sample_name}_genomad/${sample_name}_plasmid_summary.tsv", optional: true, emit: genomad_plasmid_summary
      path "${sample_name}/${sample_name}_genomad/${sample_name}_virus.fna", optional: true, emit: genomad_virus_fna
      path "${sample_name}/${sample_name}_genomad/${sample_name}_virus_summary.tsv", optional: true, emit: genomad_virus_summary
      path "${sample_name}/${sample_name}_genomad/genomad_status.txt", emit: genomad_status

    script:
    """
    set -euo pipefail
    echo "GENOMAD: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_genomad"
    TMP_PARENT=\$(mktemp -d "${sample_name}_genomad_parent.XXXXXX")
    TEMP_OUT="\$TMP_PARENT/genomad_results"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/genomad_tmp.info"
    echo "TEMP_OUT=\$TEMP_OUT" >> "\$TMP_PARENT/genomad_tmp.info"

    echo "assembly: ${assembly}" > "\$TMP_PARENT/metadata.txt"
    echo "sample: ${sample_name}" >> "\$TMP_PARENT/metadata.txt"
    echo "genomad_db: ${params.genomad_db}" >> "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.genomad_container}" >> "\$TMP_PARENT/metadata.txt"
    echo "threads: ${task.cpus}" >> "\$TMP_PARENT/metadata.txt"

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "${assembly}" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    echo "seq_bases=\$seq_bases" >> "\$TMP_PARENT/metadata.txt"

    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/genomad_status.txt"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/genomad_status.txt"
      exit 0
    fi

    genomad end-to-end -t ${task.cpus} "${assembly}" "\$TEMP_OUT" "${params.genomad_db}" > "\$TMP_PARENT/genomad.stdout" 2> "\$TMP_PARENT/genomad.stderr" || rc=\$?

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      mkdir -p "\$TARGET_DIR"
      echo "GENOMAD_FAILED rc=\$rc" > "\$TARGET_DIR/genomad_status.txt"
      echo "See tmp stderr at: \$TMP_PARENT/genomad.stderr" >> "\$TARGET_DIR/genomad_status.txt"
      cp -a "\$TMP_PARENT/genomad.stdout" "\$TMP_PARENT/genomad.stderr" "\$TARGET_DIR/" || true
      if [ -d "\$TEMP_OUT" ]; then
        cp -a "\$TEMP_OUT/." "\$TARGET_DIR/" || true
      fi
      exit 0
    fi

    # Now collect only the desired summary files into the genomad main folder
    mkdir -p "\$TARGET_DIR"

    if [ -d "\$TEMP_OUT" ]; then
      # patterns to collect
      patterns=(
        "${sample_name}_plasmid.fna"
        "${sample_name}_plasmid_summary.tsv"
        "${sample_name}_virus.fna"
        "${sample_name}_virus_summary.tsv"
      )

      for p in "\${patterns[@]}"; do
        found=\$(find "\$TEMP_OUT" -type f -name "\$p" -print -quit || true)
        if [ -n "\$found" ]; then
          cp -a "\$found" "\$TARGET_DIR/" || true
        fi
      done

      # save stdout/stderr for diagnostics
      cp -a "\$TMP_PARENT/genomad.stdout" "\$TARGET_DIR/" || true
      cp -a "\$TMP_PARENT/genomad.stderr" "\$TARGET_DIR/" || true

    else
      echo "GENOMAD_FAILED: output directory not created: \$TEMP_OUT" > "\$TARGET_DIR/genomad_status.txt"
    fi

    if [ ! -s "\$TARGET_DIR/genomad_status.txt" ]; then
      if [ ! -s "\$TARGET_DIR/${sample_name}_plasmid.fna" ] && [ ! -s "\$TARGET_DIR/${sample_name}_virus.fna" ]; then
        echo "GENOMAD_NO_HITS" > "\$TARGET_DIR/genomad_status.txt"
      else
        echo "GENOMAD_OK" > "\$TARGET_DIR/genomad_status.txt"
      fi
    fi
    """
}

/* ----------------------
   MOB-RECON per-sample process
   ---------------------- */
process MOB_RECON {
    tag { sample_name }
    container params.mob_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      path "${sample_name}/${sample_name}_mob/**", emit: mob_out
      path "${sample_name}/${sample_name}_mob/mob_status.txt", emit: mob_status
      path "${sample_name}/${sample_name}_mob/${sample_name}_mobtyper_results.txt", emit: mobtyper_results

    script:
    """
    set -euo pipefail
    echo "MOB_RECON: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_mob"
    TMP_PARENT=\$(mktemp -d "${sample_name}_mob_parent.XXXXXX")
    TEMP_OUT="\$TMP_PARENT/mob_results"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/mob_tmp.info"
    echo "TEMP_OUT=\$TEMP_OUT" >> "\$TMP_PARENT/mob_tmp.info"

    echo "assembly: ${assembly}" > "\$TMP_PARENT/metadata.txt"
    echo "sample: ${sample_name}" >> "\$TMP_PARENT/metadata.txt"
    echo "mob_db: ${params.mob_db}" >> "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.mob_container}" >> "\$TMP_PARENT/metadata.txt"
    echo "threads: ${task.cpus}" >> "\$TMP_PARENT/metadata.txt"

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "${assembly}" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    echo "seq_bases=\$seq_bases" >> "\$TMP_PARENT/metadata.txt"

    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/mob_status.txt"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/mob_status.txt"
      exit 0
    fi

    mob_recon --infile "${assembly}" --outdir "\$TEMP_OUT" --num_threads ${task.cpus} --database_directory "${params.mob_db}" > "\$TMP_PARENT/mob.stdout" 2> "\$TMP_PARENT/mob.stderr" || rc=\$?

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      mkdir -p "\$TARGET_DIR"
      echo "MOB_FAILED rc=\$rc" > "\$TARGET_DIR/mob_status.txt"
      echo "See tmp stderr at: \$TMP_PARENT/mob.stderr" >> "\$TARGET_DIR/mob_status.txt"
      if [ -d "\$TEMP_OUT" ]; then
        cp -a "\$TEMP_OUT/." "\$TARGET_DIR/" || true
      fi
      : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
      cp -a "\$TMP_PARENT/mob.stdout" "\$TARGET_DIR/" || true
      cp -a "\$TMP_PARENT/mob.stderr" "\$TARGET_DIR/" || true
      exit 0
    fi

    if [ -d "\$TEMP_OUT" ]; then
      rm -rf "\$TARGET_DIR" || true
      mv "\$TEMP_OUT" "\$TARGET_DIR"

      if [ -f "\$TARGET_DIR/chromosome.fasta" ]; then
        rm -f "\$TARGET_DIR/chromosome.fasta" || true
      fi

      if [ -f "\$TARGET_DIR/mobtyper_results.txt" ]; then
        mv "\$TARGET_DIR/mobtyper_results.txt" "\$TARGET_DIR/${sample_name}_mobtyper_results.txt" || true
        echo "MOB_OK" > "\$TARGET_DIR/mob_status.txt"
      else
        : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
        echo "MOB_NO_TYPER_RESULT" > "\$TARGET_DIR/mob_status.txt"
      fi

    else
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
      echo "MOB_FAILED: mob_recon did not create expected output dir \$TEMP_OUT" > "\$TARGET_DIR/mob_status.txt"
    fi

    [ -f "\$TARGET_DIR/${sample_name}_mobtyper_results.txt" ] || : > "\$TARGET_DIR/${sample_name}_mobtyper_results.txt"
    cp -a "\$TMP_PARENT/mob.stdout" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/mob.stderr" "\$TARGET_DIR/" || true
    """
}

/* ----------------------
   STARAMR_RUN per-sample process (keep only summary.tsv renamed)
   ---------------------- */
process STARAMR_RUN {
    tag { sample_name }
    container params.staramr_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      path "${sample_name}/${sample_name}_staramr/${sample_name}_summary.tsv", emit: staramr_summary
      path "${sample_name}/${sample_name}_staramr/staramr_status.txt", emit: staramr_status

    script:
    """
    set -euo pipefail
    echo "STARAMR_RUN: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_staramr"
    TMP_PARENT=\$(mktemp -d "${sample_name}_staramr_parent.XXXXXX")
    TMP_OUT="\$TMP_PARENT/staramr_out_\$(date +%s%N)"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/staramr_tmp.info"
    echo "TMP_OUT=\$TMP_OUT" >> "\$TMP_PARENT/staramr_tmp.info"
    echo "TARGET_DIR=\$TARGET_DIR" >> "\$TMP_PARENT/staramr_tmp.info"

    IN_ASM="${assembly}"
    echo "staramr_db: ${params.staramr_db}" > "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.staramr_container}" >> "\$TMP_PARENT/metadata.txt"

    if [ ! -s "\$IN_ASM" ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_summary.tsv"
      echo "INVALID_EMPTY_ASSEMBLY_OR_MISSING_INPUT: \$IN_ASM" > "\$TARGET_DIR/staramr_status.txt"
      exit 0
    fi

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "\$IN_ASM" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_summary.tsv"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/staramr_status.txt"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_summary.tsv"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/staramr_status.txt"
      exit 0
    fi

    if [ -e "\$TMP_OUT" ]; then
      TMP_OUT="\$TMP_PARENT/staramr_out_\$(date +%s%N)_\$RANDOM"
      echo "Adjusted TMP_OUT=\$TMP_OUT" >> "\$TMP_PARENT/staramr_tmp.info"
    fi

    echo "Running staramr search -> \$TMP_OUT" >> "\$TMP_PARENT/staramr_tmp.info"

    staramr search --database "${params.staramr_db}" --nprocs ${task.cpus} -o "\$TMP_OUT" "\$IN_ASM" > "\$TMP_PARENT/staramr.stdout" 2> "\$TMP_PARENT/staramr.stderr" || rc=\$?

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_summary.tsv"
      cp -a "\$TMP_PARENT/staramr.stdout" "\$TARGET_DIR/" || true
      cp -a "\$TMP_PARENT/staramr.stderr" "\$TARGET_DIR/" || true
      cp -a "\$TMP_PARENT/metadata.txt" "\$TARGET_DIR/" || true

      if grep -qi 'already exists' "\$TMP_PARENT/staramr.stderr" 2>/dev/null; then
        echo "STARAMR_FAILED_ALREADY_EXISTS rc=\$rc" > "\$TARGET_DIR/staramr_status.txt"
        echo "staramr stderr contains 'already exists' — TMP_OUT may have existed unexpectedly" >> "\$TARGET_DIR/staramr_status.txt"
      else
        echo "STARAMR_FAILED rc=\$rc" > "\$TARGET_DIR/staramr_status.txt"
      fi
      exit 0
    fi

    mkdir -p "\$TARGET_DIR"

    if [ -d "\$TMP_OUT" ]; then
      # copy only summary.tsv and rename
      if [ -f "\$TMP_OUT/summary.tsv" ]; then
        cp -a "\$TMP_OUT/summary.tsv" "\$TARGET_DIR/${sample_name}_summary.tsv" || true
      else
        # try to find summary.tsv anywhere below TMP_OUT
        found=\$(find "\$TMP_OUT" -type f -name 'summary.tsv' -print -quit || true)
        if [ -n "\$found" ]; then
          cp -a "\$found" "\$TARGET_DIR/${sample_name}_summary.tsv" || true
        else
          : > "\$TARGET_DIR/${sample_name}_summary.tsv"
          echo "STARAMR_FAILED: summary.tsv not found" > "\$TARGET_DIR/staramr_status.txt"
          cp -a "\$TMP_PARENT/staramr.stdout" "\$TARGET_DIR/" || true
          cp -a "\$TMP_PARENT/staramr.stderr" "\$TARGET_DIR/" || true
          exit 0
        fi
      fi
    else
      : > "\$TARGET_DIR/${sample_name}_summary.tsv"
      echo "STARAMR_FAILED: expected output directory not created: \$TMP_OUT" > "\$TARGET_DIR/staramr_status.txt"
      cp -a "\$TMP_PARENT/staramr.stdout" "\$TARGET_DIR/" || true
      cp -a "\$TMP_PARENT/staramr.stderr" "\$TARGET_DIR/" || true
      exit 0
    fi

    cp -a "\$TMP_PARENT/staramr.stdout" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/staramr.stderr" "\$TARGET_DIR/" || true

    echo "STARAMR_OK" > "\$TARGET_DIR/staramr_status.txt"
    """
}

/* ----------------------
   RGI per-sample process
   ---------------------- */
process RGI_RUN {
    tag { sample_name }
    container params.rgi_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      path "${sample_name}/${sample_name}_rgi/**", emit: rgi_out
      path "${sample_name}/${sample_name}_rgi/rgi_status.txt", emit: rgi_status
      path "${sample_name}/${sample_name}_rgi/${sample_name}_rgi_output.txt", emit: rgi_txt

    script:
    """
    set -euo pipefail
    echo "RGI_RUN: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_rgi"
    TMP_PARENT=\$(mktemp -d "${sample_name}_rgi_parent.XXXXXX")
    TMP_OUT="\$TMP_PARENT/rgi_out"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/rgi_tmp.info"
    echo "TMP_OUT=\$TMP_OUT" >> "\$TMP_PARENT/rgi_tmp.info"

    echo "assembly: ${assembly}" > "\$TMP_PARENT/metadata.txt"
    echo "sample: ${sample_name}" >> "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.rgi_container}" >> "\$TMP_PARENT/metadata.txt"
    echo "threads: ${task.cpus}" >> "\$TMP_PARENT/metadata.txt"

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "${assembly}" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    echo "seq_bases=\$seq_bases" >> "\$TMP_PARENT/metadata.txt"

    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_rgi_output.txt"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/rgi_status.txt"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_rgi_output.txt"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/rgi_status.txt"
      exit 0
    fi

    mkdir -p "\$TMP_OUT"

    rgi main -i "${assembly}" -o "\$TMP_OUT/${sample_name}_rgi_output" -n ${task.cpus} > "\$TMP_PARENT/rgi.stdout" 2> "\$TMP_PARENT/rgi.stderr" || rc=\$?

    mkdir -p "\$TARGET_DIR"
    cp -a "\$TMP_PARENT/rgi.stdout" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/rgi.stderr" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/rgi_tmp.info" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/metadata.txt" "\$TARGET_DIR/" || true

    # copy any matching outputs (including .txt)
    found=\$(find "\$TMP_OUT" -type f -name "${sample_name}_rgi_output*" -print -quit || true)
    if [ -n "\$found" ]; then
      while IFS= read -r f; do
        cp -a "\$f" "\$TARGET_DIR/" || true
      done < <(find "\$TMP_OUT" -type f -name "${sample_name}_rgi_output*" -print)
    else
      : > "\$TARGET_DIR/${sample_name}_rgi_output.txt"
      if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
        echo "RGI_FAILED rc=\$rc; output not found" > "\$TARGET_DIR/rgi_status.txt"
      else
        echo "RGI_FAILED: output not found" > "\$TARGET_DIR/rgi_status.txt"
      fi
      exit 0
    fi

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      echo "RGI_FAILED rc=\$rc" > "\$TARGET_DIR/rgi_status.txt"
    elif [ ! -s "\$TARGET_DIR/rgi_status.txt" ]; then
      echo "RGI_OK" > "\$TARGET_DIR/rgi_status.txt"
    fi
    """
}

/* ----------------------
   AMRFINDER per-sample process
   ---------------------- */
process AMRFINDER_RUN {
    tag { sample_name }
    container params.amrfinder_container

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus params.cpus
    memory params.memory
    time params.time

    input:
      tuple val(sample_name), path(assembly)

    output:
      path "${sample_name}/${sample_name}_amrfinder/**", emit: amrfinder_out
      path "${sample_name}/${sample_name}_amrfinder/${sample_name}_amrfinder_plus", emit: amrfinder_result
      path "${sample_name}/${sample_name}_amrfinder/amrfinder_status.txt", emit: amrfinder_status

    script:
    """
    set -euo pipefail
    echo "AMRFINDER_RUN: sample=${sample_name}"

    PARENT_DIR="${sample_name}"
    TARGET_DIR="\${PARENT_DIR}/${sample_name}_amrfinder"
    TMP_PARENT=\$(mktemp -d "${sample_name}_amrfinder_parent.XXXXXX")
    TMP_OUT="\$TMP_PARENT/amrfinder_out"

    mkdir -p "\$PARENT_DIR"
    echo "TMP_PARENT=\$TMP_PARENT" > "\$TMP_PARENT/amrfinder_tmp.info"
    echo "TMP_OUT=\$TMP_OUT" >> "\$TMP_PARENT/amrfinder_tmp.info"

    echo "assembly: ${assembly}" > "\$TMP_PARENT/metadata.txt"
    echo "sample: ${sample_name}" >> "\$TMP_PARENT/metadata.txt"
    echo "amrfinder_db: ${params.amrfinder_db}" >> "\$TMP_PARENT/metadata.txt"
    echo "container: ${params.amrfinder_container}" >> "\$TMP_PARENT/metadata.txt"
    echo "threads: ${task.cpus}" >> "\$TMP_PARENT/metadata.txt"

    seq_bases=\$(awk 'BEGIN{nb=0} /^>/ {next} {nb+=length(\$0)} END{print nb}' "${assembly}" || true)
    if [ -z "\$seq_bases" ]; then seq_bases=0; fi
    echo "seq_bases=\$seq_bases" >> "\$TMP_PARENT/metadata.txt"

    if [ "\$seq_bases" -eq 0 ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_amrfinder_plus"
      echo "INVALID_EMPTY_ASSEMBLY" > "\$TARGET_DIR/amrfinder_status.txt"
      exit 0
    fi
    if [ ${params.min_assembly_bases} -gt 0 ] && [ "\$seq_bases" -lt ${params.min_assembly_bases} ]; then
      mkdir -p "\$TARGET_DIR"
      : > "\$TARGET_DIR/${sample_name}_amrfinder_plus"
      echo "SKIPPED_MIN_ASSEMBLY_BASES: seq_bases=\$seq_bases min_assembly_bases=${params.min_assembly_bases}" > "\$TARGET_DIR/amrfinder_status.txt"
      exit 0
    fi

    mkdir -p "\$TMP_OUT"

    amrfinder -n "${assembly}" -o "\$TMP_OUT/${sample_name}_amrfinder_plus" --plus --threads ${task.cpus} --database "${params.amrfinder_db}" > "\$TMP_PARENT/amrfinder.stdout" 2> "\$TMP_PARENT/amrfinder.stderr" || rc=\$?

    mkdir -p "\$TARGET_DIR"
    cp -a "\$TMP_PARENT/amrfinder.stdout" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/amrfinder.stderr" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/amrfinder_tmp.info" "\$TARGET_DIR/" || true
    cp -a "\$TMP_PARENT/metadata.txt" "\$TARGET_DIR/" || true

    if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
      echo "AMRFINDER_FAILED rc=\$rc" > "\$TARGET_DIR/amrfinder_status.txt"
    fi

    if [ -f "\$TMP_OUT/${sample_name}_amrfinder_plus" ]; then
      mv "\$TMP_OUT/${sample_name}_amrfinder_plus" "\$TARGET_DIR/${sample_name}_amrfinder_plus" || true
    elif [ -f "\$TMP_OUT/${sample_name}_amrfinder_plus.tsv" ]; then
      mv "\$TMP_OUT/${sample_name}_amrfinder_plus.tsv" "\$TARGET_DIR/${sample_name}_amrfinder_plus" || true
    else
      found=\$(find "\$TMP_OUT" -type f -name "${sample_name}_amrfinder_plus*" -print -quit || true)
      if [ -n "\$found" ]; then
        cp -a "\$found" "\$TARGET_DIR/${sample_name}_amrfinder_plus" || true
      else
        : > "\$TARGET_DIR/${sample_name}_amrfinder_plus"
        if [ -n "\${rc:-}" ] && [ "\$rc" -ne 0 ]; then
          echo "AMRFINDER_FAILED rc=\$rc; output not found" > "\$TARGET_DIR/amrfinder_status.txt"
        else
          echo "AMRFINDER_FAILED: output not found" > "\$TARGET_DIR/amrfinder_status.txt"
        fi
        exit 0
      fi
    fi

    if [ ! -s "\$TARGET_DIR/amrfinder_status.txt" ]; then
      echo "AMRFINDER_OK" > "\$TARGET_DIR/amrfinder_status.txt"
    fi
    """
}

/* ----------------------
   Aggregation: collate STARAMR summaries across samples
   ---------------------- */
process COLLATE_STARAMR {
    tag "collate_staramr"

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus 1
    memory '1 GB'
    time '10m'

    input:
      path(staramr_files)

    output:
      path "${params.run_name}_staramr_summary.tsv", emit: run_staramr_summary

    script:
    """
    set -euo pipefail
    OUT="${params.run_name}_staramr_summary.tsv"

    files=( ${staramr_files.join(' ')} )
    non_empty_files=()

    for f in "\${files[@]}"; do
      if [ -s "\$f" ]; then
        non_empty_files+=( "\$f" )
      fi
    done

    if [ \${#non_empty_files[@]} -eq 0 ]; then
      echo "" > "\$OUT"
      exit 0
    fi

    head -n1 "\${non_empty_files[0]}" > "\$OUT"
    for f in "\${non_empty_files[@]}"; do
      tail -n +2 "\$f" >> "\$OUT" || true
    done
    """
}

/* ----------------------
   Aggregation: collate MOB mobtyper results across samples
   ---------------------- */
process COLLATE_MOBTYper {
    tag "collate_mobtyper"

    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus 1
    memory '1 GB'
    time '10m'

    input:
      path(mob_files)

    output:
      path "${params.run_name}_mobtyper_results.tsv", emit: run_mobtyper_summary

    script:
    """
    set -euo pipefail
    OUT="${params.run_name}_mobtyper_results.tsv"

    files=( ${mob_files.join(' ')} )
    non_empty_files=()

    for f in "\${files[@]}"; do
      if [ -s "\$f" ]; then
        non_empty_files+=( "\$f" )
      fi
    done

    if [ \${#non_empty_files[@]} -eq 0 ]; then
      echo "" > "\$OUT"
      exit 0
    fi

    head -n1 "\${non_empty_files[0]}" > "\$OUT"
    for f in "\${non_empty_files[@]}"; do
      tail -n +2 "\$f" >> "\$OUT" || true
    done
    """
}


/* ----------------------
   Aggregation: ARG summary across samples (AMRFINDER + RGI(.txt) + STARAMR)
   Option 2: Embedded format — gene(identity,coverage)
   - RGI:        Best_Hit_ARO, Best_Identities,
                 Percentage Length of Reference Sequence
                 format: gene(identity,coverage)
   - AMRFinder:  Element symbol (only where Type == "AMR"),
                 % Identity to reference, % Coverage of reference
                 format: gene(identity,coverage)
   - STARAMR:    col 3  = Genotype  (no metrics available; kept as-is)
   ---------------------- */
process COLLATE_ARG_SUMMARY {
    tag "collate_arg_summary"
    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true
    cpus 1
    memory '2 GB'
    time '20m'

    input:
      path(staramr_files)
      path(amrfinder_files)
      path(rgi_files)

    output:
      path "${params.run_name}_ARG_summary.tsv", emit: run_arg_summary

    script:
    """
    set -euo pipefail
    OUT="${params.run_name}_ARG_summary.tsv"

    echo "=== COLLATE_ARG_SUMMARY DEBUG ===" >&2
    echo "Working directory: \$(pwd)" >&2
    echo "" >&2

    echo "STARAMR files:" >&2
    for f in ${staramr_files.join(' ')}; do
        echo "  \$f" >&2
        if [[ -f "\$f" ]]; then
            echo "    Size: \$(wc -l < "\$f") lines" >&2
        fi
    done
    echo "" >&2

    echo "AMRFinder files:" >&2
    for f in ${amrfinder_files.join(' ')}; do
        echo "  \$f" >&2
        if [[ -f "\$f" ]]; then
            echo "    Size: \$(wc -l < "\$f") lines" >&2
        fi
    done
    echo "" >&2

    echo "RGI files:" >&2
    for f in ${rgi_files.join(' ')}; do
        echo "  \$f" >&2
        if [[ -f "\$f" ]]; then
            echo "    Size: \$(wc -l < "\$f") lines" >&2
        fi
    done
    echo "" >&2

    echo -e "Genome\\tRGI_ARGs\\tAMRFinder_ARGs\\tSTARAMR_ARGs" > "\$OUT"

    # Build list of samples from STARAMR files
    declare -A samples
    for f in ${staramr_files.join(' ')}; do
        if [[ -f "\$f" ]]; then
            fname=\$(basename "\$f")
            sample=\${fname%%_summary.tsv}
            if [[ -n "\$sample" ]] && [[ "\$sample" != "." ]]; then
                samples["\$sample"]=1
                echo "Found sample from STARAMR: \$sample" >&2
            fi
        fi
    done

    echo "" >&2
    echo "Total samples to process: \${#samples[@]}" >&2
    echo "" >&2

    for sample in "\${!samples[@]}"; do
        echo "--- Processing: \$sample ---" >&2
        rgi_args=""
        amr_args=""
        star_args=""

        # -------------------------------------------------------------------
        # RGI — select fields by header so column order changes do not corrupt values.
        # Output format: gene(identity,coverage)
        # -------------------------------------------------------------------
        rgi_found=false
        for f in ${rgi_files.join(' ')}; do
            if [[ -f "\$f" ]]; then
                fname=\$(basename "\$f")
                if [[ "\$fname" == "\${sample}_rgi_output.txt" ]]; then
                    echo "  RGI file: \$f" >&2
                    rgi_found=true
                    rgi_args=\$(awk -F'\\t' '
                        NR==1 {
                            for (i=1; i<=NF; i++) col[\$i]=i
                            if (!col["Best_Hit_ARO"] || !col["Best_Identities"] || !col["Percentage Length of Reference Sequence"]) {
                                print "RGI output is missing expected columns" > "/dev/stderr"
                                exit 2
                            }
                            next
                        }
                        {
                            gene = \$(col["Best_Hit_ARO"])
                            identity = \$(col["Best_Identities"])
                            coverage = \$(col["Percentage Length of Reference Sequence"])
                        }
                        gene!="" {
                            printf "%s(%s,%s)\\n", gene, identity, coverage
                        }
                    ' "\$f" | paste -sd ';')
                    gene_count=\$(echo "\$rgi_args" | awk -F';' '{print NF}')
                    echo "    Extracted \$gene_count RGI genes (with identity,coverage)" >&2
                    break
                fi
            fi
        done
        [[ "\$rgi_found" == "false" ]] && echo "  WARNING: RGI file not found for \$sample" >&2

        # -------------------------------------------------------------------
        # AMRFinder — select fields by header; retain only Type=AMR rows.
        # Output format: gene(identity,coverage)
        # -------------------------------------------------------------------
        amr_found=false
        for f in ${amrfinder_files.join(' ')}; do
            if [[ -f "\$f" ]]; then
                fname=\$(basename "\$f")
                if [[ "\$fname" == "\${sample}_amrfinder_plus" ]]; then
                    echo "  AMRFinder file: \$f" >&2
                    amr_found=true
                    amr_args=\$(awk -F'\\t' '
                        NR==1 {
                            for (i=1; i<=NF; i++) col[\$i]=i
                            coverage_col = col["% Coverage of reference"] ? col["% Coverage of reference"] : col["% Coverage of reference sequence"]
                            identity_col = col["% Identity to reference"] ? col["% Identity to reference"] : col["% Identity to reference sequence"]
                            if (!col["Element symbol"] || !col["Type"] || !coverage_col || !identity_col) {
                                print "AMRFinder output is missing expected columns" > "/dev/stderr"
                                exit 2
                            }
                            next
                        }
                        \$(col["Type"])=="AMR" && \$(col["Element symbol"])!="" {
                            gene     = \$(col["Element symbol"])
                            coverage = \$(coverage_col)
                            identity = \$(identity_col)
                            printf "%s(%s,%s)\\n", gene, identity, coverage
                        }
                    ' "\$f" | paste -sd ';')
                    gene_count=\$(echo "\$amr_args" | awk -F';' '{print NF}')
                    echo "    Extracted \$gene_count AMR genes (Type=AMR, with identity,coverage)" >&2

                    # Show Type column distribution for debugging
                    type_counts=\$(awk -F'\\t' 'NR==1 {for (i=1; i<=NF; i++) col[\$i]=i; next} col["Type"] {print \$(col["Type"])}' "\$f" | sort | uniq -c)
                    echo "    Type distribution:" >&2
                    echo "\$type_counts" | sed 's/^/      /' >&2
                    break
                fi
            fi
        done
        [[ "\$amr_found" == "false" ]] && echo "  WARNING: AMRFinder file not found for \$sample" >&2

        # -------------------------------------------------------------------
        # STARAMR — select Genotype by header; split comma-separated genes, join with ;
        # No metrics available from summary.tsv — kept as plain gene names
        # -------------------------------------------------------------------
        star_found=false
        for f in ${staramr_files.join(' ')}; do
            if [[ -f "\$f" ]]; then
                fname=\$(basename "\$f")
                if [[ "\$fname" == "\${sample}_summary.tsv" ]]; then
                    echo "  STARAMR file: \$f" >&2
                    star_found=true
                    star_args=\$(awk -F'\\t' '
                        NR==1 {
                            for (i=1; i<=NF; i++) col[\$i]=i
                            if (!col["Genotype"]) {
                                print "StarAMR output is missing Genotype column" > "/dev/stderr"
                                exit 2
                            }
                            next
                        }
                        \$(col["Genotype"])!="" {print \$(col["Genotype"])}
                    ' "\$f" | tr ',' ';' | tr -d ' ')
                    gene_count=\$(echo "\$star_args" | awk -F';' '{print NF}')
                    echo "    Extracted \$gene_count STARAMR genes" >&2
                    break
                fi
            fi
        done
        [[ "\$star_found" == "false" ]] && echo "  WARNING: STARAMR file not found for \$sample" >&2

        echo -e "\${sample}\\t\${rgi_args}\\t\${amr_args}\\t\${star_args}" >> "\$OUT"
        echo "" >&2
    done

    echo "=== Final Output ===" >&2
    cat "\$OUT" >&2
    echo "===================" >&2
    """
}

process WRITE_TOOL_STATUS_SUMMARY {
    tag "write_tool_status_summary"
    publishDir "${params.outDir}/${params.run_name}", mode: 'copy', overwrite: true

    cpus 1
    memory '1 GB'
    time '10m'

    input:
      tuple val(sample_names), path(status_files, stageAs: 'status??/*')
      path(report_files)

    output:
      path "${params.run_name}_tool_status_summary.tsv", emit: tool_status_summary

    script:
    """
    set -euo pipefail
    OUT="${params.run_name}_tool_status_summary.tsv"
    printf "Sample\\tStatus file\\tResult\\n" > "\$OUT"
    files=( ${status_files.join(' ')} )
    samples=( ${sample_names.join(' ')} )
    for i in "\${!files[@]}"; do
      result=\$(head -n1 "\${files[\$i]}")
      if [ -z "\$result" ]; then
        echo "Empty status file: \${files[\$i]}" >&2
        exit 1
      fi
      printf "%s\\t%s\\t%s\\n" "\${samples[\$i]}" "\$(basename "\${files[\$i]}")" "\$result"
    done | LC_ALL=C sort >> "\$OUT"
    """
}

process FAIL_ON_TOOL_ERRORS {
    tag "fail_on_tool_errors"

    cpus 1
    memory '1 GB'
    time '10m'

    input:
      path(status_summary)

    script:
    """
    set -euo pipefail
    if awk -F '\\t' 'NR > 1 && (\$3 ~ /FAILED|INVALID/ || \$3 ~ /^(WARNING:|ERROR:)/) { found=1 } END { exit !found }' "${status_summary}"; then
      echo "One or more tools failed; see the published status summary and per-sample diagnostics." >&2
      exit 1
    fi
    """
}

/* ----------------------
   Workflow entrypoint
   ---------------------- */
workflow {
    def bakta_results    = BAKTA_ANNOTATE(assemblies_ch)

    def genomad_results  = GENOMAD_RUN(assemblies_ch)
    def mob_results      = MOB_RECON(assemblies_ch)
    def staramr_results  = STARAMR_RUN(assemblies_ch)

    // New tools
    def rgi_results      = RGI_RUN(assemblies_ch)
    def amrfinder_results= AMRFINDER_RUN(assemblies_ch)

    // Collate per-tool
    def staramr_collated = COLLATE_STARAMR(staramr_results.staramr_summary.collect())
    def mob_collated = COLLATE_MOBTYper(mob_results.mobtyper_results.collect())

    // ARG summary table (TSV format, Option 2: embedded metrics)
    def arg_collated = COLLATE_ARG_SUMMARY(
      staramr_results.staramr_summary.collect(),
      amrfinder_results.amrfinder_result.collect(),
      rgi_results.rgi_txt.collect()
    )

    def all_tool_statuses = bakta_results.bakta_status
        .mix(genomad_results.genomad_status)
        .mix(mob_results.mob_status)
        .mix(staramr_results.staramr_status)
        .mix(rgi_results.rgi_status)
        .mix(amrfinder_results.amrfinder_status)
        .collect()
    def status_inputs = all_tool_statuses.map { files ->
        tuple(files.collect { it.parent.parent.name }, files)
    }
    def all_reports = staramr_collated.run_staramr_summary
        .mix(mob_collated.run_mobtyper_summary)
        .mix(arg_collated.run_arg_summary)
        .collect()
    def tool_status_summary = WRITE_TOOL_STATUS_SUMMARY(status_inputs, all_reports)
    FAIL_ON_TOOL_ERRORS(tool_status_summary.tool_status_summary)
}

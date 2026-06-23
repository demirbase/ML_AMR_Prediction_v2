#!/usr/bin/env nextflow
// ============================================================================
// 08_blast_pipeline.nf  —  BLAST Annotation Pipeline
// ============================================================================
// Nextflow pipeline that takes the top-k-mer FASTA generated in Step 07 and
// runs two parallel annotation passes:
//   1. CARD local BLAST   — queries the Comprehensive Antibiotic Resistance
//                           Database to identify known resistance gene k-mers.
//   2. NCBI remote BLAST  — queries the full NCBI nt database remotely for
//                           any uncharacterised hits.
//
// Both processes produce tabular (outfmt 6) TSV output published into the
// 05_explainability subdirectory.
//
// Usage (invoked by 08_blast_annotation.py):
//   nextflow run scripts/08_blast_pipeline.nf \
//     --fasta        path/to/02_top_features_{antibiotic}.fasta \
//     --card_db      path/to/data/blast_db/card_nt/card \
//     --outdir       path/to/analysis_results/{antibiotic}/05_explainability \
//     --antibiotic   ciprofloxacin \
//     --threads      8 \
//     --evalue       10 \
//     --word_size    11
// ============================================================================

nextflow.enable.dsl = 2

// ---------------------------------------------------------------------------
// Parameters (overridable from CLI or 08_blast_annotation.py subprocess call)
// ---------------------------------------------------------------------------
params.fasta      = ""
params.card_db    = ""
params.antibiotic = "unknown"
params.outdir     = "results/${params.antibiotic}/05_explainability"
params.threads    = 8
params.evalue     = 10
params.word_size  = 11
// BLAST task: 'blastn-short' for short queries (<50 bp — k-mers AND short
// unitigs), 'blastn' for longer. 08_blast_annotation.py sets this from the
// ACTUAL longest query in the FASTA. word_size 11 is sensitive for both.
params.task       = "blastn-short"

// ---------------------------------------------------------------------------
// Shared BLAST output format: standard tabular + extra annotation fields
// ---------------------------------------------------------------------------
// qlen (query length) is emitted before stitle so 09 can compute coverage as
// alignment_length / query_length — correct for variable-length unitig queries
// (k=21 is wrong for unitigs). stitle stays last (free-text, tab-safe).
// NOTE: a params assignment (not a top-level `def`) — newer Nextflow's strict
// parser rejects top-level statements mixed with process declarations.
params.outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle"

// ---------------------------------------------------------------------------
// PROCESS 1: Local BLAST against CARD database
// ---------------------------------------------------------------------------
process CARD_BLAST {
    tag          "CARD | ${params.antibiotic}"
    publishDir   params.outdir, mode: 'copy', overwrite: true
    errorStrategy 'ignore'   // Missing DB must not abort the NCBI process

    input:
    path fasta

    output:
    path "03_card_blast_results_${params.antibiotic}.tsv"

    script:
    """
    blastn \\
        -query       ${fasta} \\
        -db          ${params.card_db} \\
        -task        ${params.task} \\
        -dust        no \\
        -outfmt      "${params.outfmt}" \\
        -evalue      ${params.evalue} \\
        -word_size   ${params.word_size} \\
        -num_threads ${params.threads} \\
        -out         03_card_blast_results_${params.antibiotic}.tsv
    """
}

// ---------------------------------------------------------------------------
// PROCESS 2: Remote BLAST against NCBI nt (no local DB required)
// ---------------------------------------------------------------------------
process NCBI_REMOTE_BLAST {
    tag           "NCBI remote | ${params.antibiotic}"
    publishDir    params.outdir, mode: 'copy', overwrite: true
    errorStrategy 'ignore'   // A remote/network failure must not abort the CARD process

    input:
    path fasta

    output:
    path "04_ncbi_blast_results_${params.antibiotic}.tsv"

    script:
    """
    blastn \\
        -query     ${fasta} \\
        -db        nt \\
        -remote \\
        -task      ${params.task} \\
        -dust      no \\
        -outfmt    "${params.outfmt}" \\
        -evalue    ${params.evalue} \\
        -word_size ${params.word_size} \\
        -out       04_ncbi_blast_results_${params.antibiotic}.tsv
    """
}

// ---------------------------------------------------------------------------
// WORKFLOW: Run both processes in parallel from the same FASTA channel
// ---------------------------------------------------------------------------
workflow {
    if (!params.fasta) {
        error "ERROR: --fasta parameter is required."
    }

    // Two independent channels so one process failure never blocks the other
    fasta_card_ch  = Channel.fromPath(params.fasta, checkIfExists: true)
    fasta_ncbi_ch  = Channel.fromPath(params.fasta, checkIfExists: true)

    CARD_BLAST(fasta_card_ch)
    NCBI_REMOTE_BLAST(fasta_ncbi_ch)
}

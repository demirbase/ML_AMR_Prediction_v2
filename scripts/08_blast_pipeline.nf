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
params.organism   = "ecoli"
// organism-scoped default (audit Issue 28); 08_blast_annotation.py passes the
// fully-resolved --outdir, so this default only affects a standalone `nextflow run`.
params.outdir     = "results/${params.organism}/${params.antibiotic}/05_explainability"
params.threads    = 8
params.evalue     = 10
params.word_size  = 11
// BLAST task: 'blastn-short' for short queries (<50 bp — k-mers AND short
// unitigs), 'blastn' for longer. 08_blast_annotation.py sets this from the
// ACTUAL longest query in the FASTA. word_size 11 is sensitive for both.
// NOTE: this 'task'/'word_size' pair drives the LOCAL CARD pass only.
params.task       = "blastn-short"

// ---------------------------------------------------------------------------
// NCBI remote pass — DECOUPLED from CARD (see NCBI_REMOTE_BLAST below)
// ---------------------------------------------------------------------------
// The public NCBI BLAST server kills 'blastn-short' + word_size 7 over the full
// nt database with SIGXCPU (CPU-usage limit) — a short 7-base seed generates far
// too many extensions across nt, even when restricted to one species. So the
// remote pass MUST use 'blastn' + a larger word_size; this is still sufficient
// because we are looking for high-identity genomic-context copies (the unitigs
// occur near-identically in nt), not weak homology. entrez_query restricts the
// remote search to the study organism (set by 08_blast_annotation.py from the
// registry display_name — never hardcoded); empty = no restriction.
params.ncbi_task        = "blastn"
params.ncbi_word_size   = 11
params.max_target_seqs  = 50
params.entrez_query     = ""

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
    // Decoupled remote params (blastn/word_size 11 — see params block above).
    // entrez_query is added only when non-empty; -num_threads is intentionally
    // omitted because '-remote' runs on NCBI's servers (local threading ignored).
    def entrez_flag = params.entrez_query ? "-entrez_query '${params.entrez_query}'" : ""
    """
    blastn \\
        -query           ${fasta} \\
        -db              nt \\
        -remote \\
        -task            ${params.ncbi_task} \\
        -dust            no \\
        -outfmt          "${params.outfmt}" \\
        -evalue          ${params.evalue} \\
        -word_size       ${params.ncbi_word_size} \\
        -max_target_seqs ${params.max_target_seqs} \\
        ${entrez_flag} \\
        -out             04_ncbi_blast_results_${params.antibiotic}.tsv
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

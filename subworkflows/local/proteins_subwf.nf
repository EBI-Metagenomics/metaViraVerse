include { PHAMMSEQS                        } from '../../modules/local/phammseqs'

include { HMMER_HMMSEARCH                  } from '../../modules/nf-core/hmmer/hmmsearch'

include { AMR_ANNOTATION                   } from '../ebi-metagenomics/amr_annotation'


workflow PROTEINS_PROCESSING {

    take:
    input  // (meta, reps_faa.uncompressed, reps_gff)

    main:

    ch_versions = channel.empty()

    ch_proteins = input.map{ meta, faa, _gff -> tuple(meta, faa) }

    //
    // -------- Antimicrobial resistence detection
    //
    AMR_ANNOTATION (
        input,
        params.amrfinderplus_db,
        params.deeparg_db,
        params.deeparg_db_version,
        params.deeparg_model,
        params.deeparg_tool_version,
        params.rgi_db,
        params.skip_amrfinderplus,
        params.skip_deeparg,
        params.skip_rgi
    )

    if (params.phammseqs) {
        //
        // -------- Assort phage protein sequences into phamilies using MMseqs2
        //
        PHAMMSEQS(
           ch_proteins
        )
        ch_versions = ch_versions.mix(PHAMMSEQS.out.versions)
        // TODO: maybe run PHAMCLUST on small dataset
    }

    //
    // -------- Functional annotation
    //
    def hmm_ch = channel
        .fromPath("${params.annotation_db}/*.hmm.gz")
        .map { hmm_file -> tuple(hmm_file.baseName.split('.')[0], hmm_file) }
    hmm_ch.view()

    def hmmsearch_input = hmm_ch
        .combine(ch_proteins)
        .map { hmm_id, hmm_file, faa_id, faa_file ->
            def meta = [
                id: "${hmm_id}",
                hmm_name: hmm_file.name,
                faa_name: faa_file.name
            ]
            tuple(meta, hmm_file, faa_file, true, true, true)
        }
    HMMER_HMMSEARCH (
        hmmsearch_input
    )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions)

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}
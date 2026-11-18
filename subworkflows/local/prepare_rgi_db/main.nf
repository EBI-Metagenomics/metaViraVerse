nextflow.enable.types = true

include { RGI_DOWNLOADDB     } from '../../../modules/ebi-metagenomics/rgi/downloaddb/main'
include { RGI_CARDANNOTATION } from '../../../modules/nf-core/rgi/cardannotation/main'

/*
 * The RGI (CARD) database, processed for RGI once per run and shared by the viral and
 * the plasmid AMR annotation:
 *   - no --rgi_db: download the latest CARD and process it with `rgi card_annotation`
 *   - --rgi_db that is already processed (path contains 'card_database_processed'): used as it is
 *   - any other --rgi_db (a directory with card.json): processed with `rgi card_annotation`
 * Null when RGI is skipped for both viruses and plasmids.
*/
workflow PREPARE_RGI_DB {

    take:
    rgi_db: String?
    skip: Boolean       // RGI skipped for viruses and for plasmids

    main:
    def processed: Value<Path?> = channel.value(null)
    if ( !skip ) {
        if ( !rgi_db ) {
            def card = RGI_DOWNLOADDB()
            processed = RGI_CARDANNOTATION(card.card_json).db
        } else if ( !rgi_db.contains('card_database_processed') ) {
            processed = RGI_CARDANNOTATION(file(rgi_db, checkIfExists: true)).db
        } else {
            processed = channel.value(file(rgi_db, checkIfExists: true))
        }
    }

    emit:
    db: Value<Path?> = processed
}

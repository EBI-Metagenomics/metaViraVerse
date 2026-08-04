process REPSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

   input:
    tuple val(meta), path(proteins)
    path(db_search)
    path(db_filter)

   output:
   tuple val(meta), path("prots_vs_rep.tsv"), emit: prots_vs_rep
   tuple val(meta), path("rep_domains.tsv"),  emit: rep_domains

   script:
   """
   echo "hmsearch"
   hmmsearch --cut_ga --cpu ${task.cpus} -o log --domtblout prots_vs_rep.tsv ${db_search} ${proteins}

   echo "DomainArch"
   plasquid_dom_arch.R prots_vs_rep.tsv  # output: multi_dom_rip.tsv, single_dom_rip.tsv, domain_architecture.RDS

   echo "FilterDomain"
   plasquid_filter_rip.R multi_dom_rip.tsv single_dom_rip.tsv domain_architecture.RDS ${db_filter}  # output: rep_domains.tsv
   """
}

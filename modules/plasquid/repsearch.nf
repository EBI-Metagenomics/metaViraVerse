process REPSEARCH {

    label 'process_medium'
    tag "${meta.id}"
    container "docker://mgimenez720/plasquid:latest"

   input:
    tuple val(meta), path(proteins)
    path(db_search)
    path(db_filter)

   output:
   tuple val(meta), path("ProtsvsRep.tsv"), emit: prots_vs_rep
   tuple val(meta), path("Rep_domains.tsv"), emit: rep_domains

   script:
   """
   echo "hmsearch"
   hmmsearch --cut_ga --cpu ${task.cpus} -o log --domtblout ProtsvsRep.tsv ${db_search} ${proteins}

   echo "DomainArch"
   Dom_Arch.R ProtsvsRep.tsv  # output: multi_dom_RIP.tsv, single_dom_RIP.tsv, Domain_Architecture.RDS

   echo "FilterDomain"
   Filter_RIP.R multi_dom_RIP.tsv single_dom_RIP.tsv Domain_Architecture.RDS ${db_filter}
   """
}

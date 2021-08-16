#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process HmPOS {

    input:
      val pgs_id

    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmPOS $pgs_id $params.genebuild -loc_files $params.loc_files -loc_hmoutput $params.loc_hmoutput --gzip
    """
}

process HmVCF {

    input:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmVCF $pgs_id $params.genebuild -loc_files $params.loc_hmoutput -loc_hmoutput $params.loc_hmoutput -loc_vcfs $params.loc_vcfs
    """
}

workflow {
    data = channel.from(params.pgs)
    HmPOS(data)
    HmVCF(HmPOS.out)
}

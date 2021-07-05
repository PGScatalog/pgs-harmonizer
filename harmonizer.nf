#!/usr/bin/env nextflow

harmonize_ch = Channel.fromList(params.pgs)

process HmPOS {

    input:
    val pgs_id from harmonize_ch

    output:
    val pgs_id into HmPOS_PGSs

    """
    python $params.loc_scripts/Harmonize.py HmPOS $pgs_id $params.genebuild -loc_files $params.loc_files -loc_hmoutput $params.loc_hmoutput --gzip
    """
}

process HmVCF {

    input:
    // val HmPOS_file from HmPOS_files
    val pgs_id from HmPOS_PGSs

    """
    python $params.loc_scripts/Harmonize.py HmVCF $pgs_id $params.genebuild -loc_files $params.loc_hmoutput -loc_hmoutput $params.loc_hmoutput -loc_vcfs $params.loc_vcfs
    """
}

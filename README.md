# pgs-harmonizer
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 


## Requirements
- Python packages: requests, pandas, pyliftover


## pseudocode (adapted from GWAS Catalog [README](https://github.com/EBISPOT/sum-stats-formatter/blob/master/harmonisation/README.md))
**Implemented**:
<pre>READ/PARSE PGS Scoring File and headers

FOR each variant in file
    IF RSID maps to genomic location in Ensembl THEN
        update locations based on Ensembl mapping
        IF RSID != original RSID THEN
            update rsID (provide original rsID in hm_info)
        check and flag if the alleles are consistent with the API
    ELIF Able to liftover locations to current build THEN
        liftover locations to current build
        check and flag if the alleles are consistent with the ENSEMBL VCF
    ELSE
        *flag* variant and provide original mappings in hm_info column as dictionary
    ENDIF
ENDFOR
</pre>

     +----+--------------------------------------------------------------+
     |Code|Description of harmonisation process                          |
     +----+--------------------------------------------------------------+
     |0   | Author-reported variant information                          |
     |1   | Mapped by rsID [allele(s) match ENSEMBL API]                 |
     |2   | Lifted over by position                                      |
     |3   | Lifted over by position [multiple mapped locations]          |
     |-1  | Unable to map the variant                                    |
     |-2  | Mapped by rsID [allele(s) map to reverse strand]             |
     |-3  | Mapped by rsID [alleles do not match ENSEMBL]                |
     |-4  | Strands flipped [reverse complement alleles exist in VCF]    |
     |-5  | Variant doesn't exist in the ENSEMBL VCF                     |
     +----+--------------------------------------------------------------+

## Related resources
### OpenTargets [GWAS summary statistics harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
- Harmonization flowchart : https://github.com/opentargets/genetics-sumstat-harmoniser/blob/master/flowchart_v3.svg
### GWAS Catalog [Summary Statistics harmonisation](https://github.com/EBISPOT/sum-stats-formatter/tree/master/harmonisation)
- Slightly different but uses the OpenTargets code/pipeline
### Liftover tools
- UCSC liftover tools 
- [pyliftover](https://pypi.org/project/pyliftover/) - code: https://github.com/konstantint/pyliftover

### ENSEMBL APIs
### Variation POST (used to find location by rsID):
- GRCh38 (Current release): https://rest.ensembl.org/documentation/info/variation_post
- GRCh37 (Past release): http://grch37.rest.ensembl.org/documentation/info/variation_post
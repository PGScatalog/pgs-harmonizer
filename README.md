# pgs-harmoniser
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 


## Requirements
- Python packages: requests, pandas


## pseudocode (adapted from GWAS Catalog)
<pre>format PGS Scoring File (e.g. headers)

FOR each variant in file
    IF RSID maps to genomic location in Ensembl THEN
        Update locations based on Ensembl mapping
    ELIF can liftover locations to current build THEN
        liftover locations to current build
    ELSE
        *flag* variant
    ENDIF
ENDFOR

FOR each non-palindromic variant
    check orientation (query Ensembl reference VCF with chr:bp, effect and other alleles)
ENDFOR

summarise the orientation of the variants: outcomes are ‘forward’, ‘reverse’ or ‘mixed’

IF ‘mixed’ THEN
    *flag* palindromic variants
ELSE
    proceed with all variants (including palindromic snps) assuming consensus orientation
ENDIF

FOR each remaining variant:
    get rsid and update variant_id
    check orientation of variant against Ensembl reference VCF
    orientate variant to reference (can flip alleles, betas, ORs, CIs allele frequencies)
    remove if variant_id, p_value, base_pair_location, chromosome are invalid
ENDFOR
</pre>
### Related resources
#### OpenTargets [GWAS summary statistics harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
- Harmonization flowchart : https://github.com/opentargets/genetics-sumstat-harmoniser/blob/master/flowchart_v3.svg
#### GWAS Catalog [Summary Statistics harmonisation](https://github.com/EBISPOT/sum-stats-formatter/tree/master/harmonisation)
- Slightly different but uses the OpenTargets code/pipeline
#### Liftover tools
- UCSC liftover tools 
- [pyliftover](https://pypi.org/project/pyliftover/) - code: https://github.com/konstantint/pyliftover

### ENSEMBL APIs
### Variation post (used to find location by rsID):
- GRCh38 (Current release): https://rest.ensembl.org/documentation/info/variation_post
- GRCh37 (Past release): http://grch37.rest.ensembl.org/documentation/info/variation_post
# pgs-harmonizer
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 


## Requirements
- Python packages: requests, pandas, pyliftover

## Current options
<pre>$ python Harmonize.py -h
usage: Harmonize.py [-h] -id PGSID -build GENOMEBUILD [-loc_scorefiles DIR]
                    [-source_build BUILD] [--var2location]
                    [--addReferenceAllele] [--ignore_rsid] [--gzip]

Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome
build.

optional arguments:
  -h, --help            show this help message and exit
  -id PGSID             PGS Catalog Score ID
  -build GENOMEBUILD    Target genome build choices = ['GRCh37', GRCh38']
  -loc_scorefiles DIR   Root directory where the PGS files are located,
                        otherwise assumed to be in: ../pgs_ScoringFiles/
  -source_build BUILD   Source genome build [overwrites the scoring file
                        header information]
  --var2location        Uses the annotations from the var2location.pl script
                        (ENSEMBL SQL connection)
  --addReferenceAllele  Adds the reference_allele(s) column for rsIDs that
                        only have a recorded effect_allele
  --ignore_rsid         Ignores rsIDs and attempts to harmonize variants by
                        liftover only
  --gzip                Writes gzipped harmonized output</pre>

## pseudocode (adapted from GWAS Catalog [README](https://github.com/EBISPOT/sum-stats-formatter/blob/master/harmonisation/README.md))
<pre>READ/PARSE PGS Scoring File and headers

FOR each variant in file
    IF RSID maps to genomic location in Ensembl THEN
        update locations based on Ensembl mapping
        IF RSID != original RSID THEN
            update rsID (provide original rsID in hm_info)
    ELIF Able to liftover locations to current build THEN
        liftover locations to current build
    ELSE
        *flag* variant and provide original mappings in hm_info column as dictionary
    ENDIF
    
    Check and flag if the alleles are consistent with the ENSEMBL VCF
    Write output
ENDFOR
</pre>

     +----+--------------------------------------------------------------+
     |Code|Description of harmonisation process                          |
     +----+--------------------------------------------------------------+
     | 5  | Variant mapped [exists in ENSEMBL VCF]                       |
     | 4  | Variant mapped [exists in ENSEMBL VCF with ambiguous         |
     |    |                  orientation (e.g. A/T, C/G variants]        |
     | 0  | Author-reported variant information                          |
     |-1  | Unable to map the variant                                    |
     |-4  | Strands flipped? [reverse complement alleles exist in VCF]   |
     |-5  | Variant doesn't exist in the ENSEMBL VCF                     |
     +----+--------------------------------------------------------------+

## Related resources
### OpenTargets [GWAS summary statistics harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
- Harmonization flowchart : https://github.com/opentargets/genetics-sumstat-harmoniser/blob/master/flowchart_v3.svg
### GWAS Catalog [Summary Statistics harmonisation](https://github.com/EBISPOT/gwas-sumstats-harmoniser)
- Slightly different but uses the OpenTargets code/pipeline
### Liftover tools
- UCSC liftover tools 
- [pyliftover](https://pypi.org/project/pyliftover/) - code: https://github.com/konstantint/pyliftover

### ENSEMBL APIs
### Variation POST (used to find location by rsID):
- GRCh38 (Current release): https://rest.ensembl.org/documentation/info/variation_post
- GRCh37 (Past release): http://grch37.rest.ensembl.org/documentation/info/variation_post
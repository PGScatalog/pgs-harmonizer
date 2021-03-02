# pgs-harmonizer
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 


## Requirements
- Python packages: requests, pandas, pyliftover, tqdm, cyvcf2

## Current options
<pre>$ python ./Harmonize.py -h
usage: Harmonize.py [-h] -id PGS###### -build GRCh## [-loc_scorefiles DIR]
                    [-source_build GENOMEBUILD] [-cohort_vcf COHORT]
                    [--var2location] [--addOtherAllele] [--ignore_rsid]
                    [--gzip]

Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome
build.

optional arguments:
  -h, --help            show this help message and exit
  -id PGS######         PGS Catalog Score ID
  -build GRCh##         Target genome build choices: 'GRCh37'or GRCh38'
  -loc_scorefiles DIR   Root directory where the PGS files are located,
                        otherwise assumed to be in: ../pgs_ScoringFiles/
  -source_build GENOMEBUILD
                        Source genome build [overwrites information in the
                        scoring file header]
  -cohort_vcf COHORT    Cohort VCF: Used to check if a variant is present in
                        the genotyped/imputed variants for a cohort and add
                        other allele when the information from ENSEMBL is
                        ambiguous (multiple potential alleles)
  --var2location        Uses the annotations from the var2location.pl script
                        (ENSEMBL SQL connection)
  --addOtherAllele      Adds a other_allele(s) column for PGS that only have a
                        recorded effect_allele
  --ignore_rsid         Ignores rsID mappings and harmonizes variants using
                        only liftover
  --skip_strandflips    This flag will stop the harmonizing from trying to
                        correct strand flips
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
    
    CHECK variant alleles against ENSEMBL or cohort-specific VCF and flag if the alleles are consistent (e.g. present, flipped, palindromic, etc)    
ENDFOR
</pre>

     +----+--------------------------------------------------------------+
     |Code|Description of harmonisation process                          |
     +----+--------------------------------------------------------------+
     | 5  | Mapped [variant exists in REFERENCE VCF]                     |
     | 4  | Mapped [variant exists in REFERENCE VCF, other allele(s) with| 
     |    |         ambiguous orientations exist at the locus]           |
     | 3  | Mapped [variant in REFERENCE VCF with ambiguous orientation  |
     |    |         (e.g. A/T, C/G variants) ]                           |
     | 0  | Author-reported variant information                          |
     |-1  | Unable to map the variant                                    |
     |-4  | Strands flipped? [reverse complement alleles exist in VCF,   |
     |    |                  default behaviour is to correct by flipping]|
     |-5  | Variant doesn't exist in the REFERENCE VCF                   |
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
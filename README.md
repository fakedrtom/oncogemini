OncoGEMINI
=============================================================================

Overview
========
OncoGEMINI is an adaptation of [GEMINI](https://github.com/arq5x/gemini) intended for the improved identification of
biologically and clincally relevant tumor variants from multi-sample and longitudinal
tumor sequencing data. Using a GEMINI-compatible database (generated from an annotated 
VCF file), OncoGEMINI is able to filter tumor variants based on included genomic
annotations and various allele frequency signatures. 


Installation
============


Documentation
================
Since OncoGEMINI retains much of the functionality of GEMINI, it may also be 
helpful to refer to GEMINI's official documentation which can be found [here](http://gemini.readthedocs.org/en/latest/).

VCF Preparation
----------------
Like GEMINI, multi-allelic sites need to be decomposed and normalized using [vt](https://genome.sph.umich.edu/wiki/Vt).
A more thorough explanation and guide for doing this can be found at the [GEMINI documentation](https://gemini.readthedocs.io/en/latest/#new-gemini-workflow).
Similarly, it is recommended that a VCF be annotated with either VEP or SnpEff before
additional annotations are included.

OncoGEMINI relies upon VCF annotations for the creation of searchable fields within 
a database. Therefore it is important that a VCF be annotated with all information that
a user desires for filtering. With that in mind, OncoGEMINI was designed to be used 
alongside [vcfanno](https://github.com/brentp/vcfanno) to accomplish all VCF annotation needs. 
Please consult the vcfanno link for details regarding its proper usage, but in short, with
a completed vcfanno configuration file, VCFs can be annotated quite simply:

```
./vcfanno vcfanno.config prepared.vcf.gz > annotated.vcf.gz
```

OncoGEMINI was also developed alongside [CRAB](https://github.com/fakedrtom/cancer_annotations) and many useful, cancer-relevant
annotations can be found and downloaded there, including a vcfanno configuration for many
of the included annotations.

Database Creation
----------------
Properly prepared and annotated VCFs can be used to create OncoGEMINI databases with [vcf2db](https://github.com/quinlan-lab/vcf2db).
The creation of a database with vcf2db also requires a pedigree-like file, referred to as a 
sample manifest, to be included. The structure of this file is similar to a more traditional
pedigree file, but inlcudes additional columns corresponding to a patient
identifier, the sequential point in which that sample was obtained (to reflect longitudinal
data across multiple timepoints where time = 0 reflect a normal or non-tumor sample and time > 0 indicates
tumor samples with different sampling times), and any sample purity values, if known.

```
#family_id      name    paternal_id     maternal_id     sex     phenotype       patient_id      time    purity
1               A0      0               0               2       1               A               0       0
1               A1      0               0               2       2               A               1       0.1
1               A2      0               0               2       2               A               2       0.3
1               B0      0               0               2       1               B               0       0
1               B1      0               0               2       2               B               1       0.5
```

Together, the annotated VCF and sample manifest file are used by the vcf2db script to generate
the OncoGEMINI database:

```
python vcf2db.py annotated.vcf.gz sample.manifest database.db
```

Usage
----------------
OncoGEMINI utilizes SQL queries in combination with tool commands to search the 
database for variants that match requested filters. 

### query
For general searches, the *query* tool allows for customization. This is carried over
from GEMINI and further details can be found [here](https://gemini.readthedocs.io/en/latest/content/querying.html#basic-queries).
For example, to search for a specific variant found on chromosome 13 at position 32,900,000, the following *query*
command would return the chromosome, start and end positions, reference and alternate alleles and gene for 
the specified variant (if it exists):

```
oncogemini query -q "select chrom, start, end, ref, alt, gene from variants where chrom == 13 and start == 32899999 and end == 32900000" database.db
```
### bottleneck
The *bottleneck* tool is designed to identify variants whose allele frequencies increase across
sampling timepoints. By default, *bottleneck* will require a variant to be absent in any normal
samples, the slope made by all included allele frequencies is greater than 0.05, and the R
correlation coefficient fo all allele frequencies if greater than 0.5, but these and other parameters
can be adjusted with the following usage options:

```
usage: gemini bottleneck [-h] [--minDP INTEGER] [--minGQ INTEGER]
                         [--maxNorm FLOAT] [--minSlope FLOAT] [--minR FLOAT]
                         [--samples STRING] [--minEnd FLOAT] [--endDiff FLOAT]
                         [--patient STRING] [--columns STRING]
                         [--filter STRING] [--purity] [--somatic_only]
                         [--cancers STRING]
                         db

positional arguments:
  db                The name of the database to be queried

optional arguments:
  -h, --help        show this help message and exit
  --minDP INTEGER   Minimum depth required in all samples default is 0)
  --minGQ INTEGER   Minimum genotype quality required in all samples (default
                    is 0)
  --maxNorm FLOAT   Specify a maximum normal sample AF to allow (default is 0)
  --minSlope FLOAT  Minimum slope required for the AFs across samples (default
                    is 0.05)
  --minR FLOAT      Minimum r correlation coefficient required for AFs
                    (default is 0.5)
  --samples STRING  Rather than including all samples, a string of comma-
                    separated specified samples to use (default is "All")
  --minEnd FLOAT    Minimum AF required of the sample representing the final
                    timepoint (default is 0)
  --endDiff FLOAT   Minimum required AF difference between the samples
                    representing the first and final timepoints (default is 0)
  --patient STRING  Specify a patient to filter (should correspond to a
                    patient_id in ped file)
  --columns STRING  A list of columns that you would like returned (default is
                    "*", which returns every column)
  --filter STRING   Restrictions to apply to variants (SQL syntax)
  --purity          Using purity estimates in ped file, make corrections to AF
                    to be used
  --somatic_only    Only include variants that have been marked as somatic via
                    the set_somatic command
  --cancers STRING  Restrict results to variants/genes associated with
                    specific cancer types by entering a comma-separated string
                    of cancer type abbreviations (see documents for
                    abbreviations) REQUIRES that db include
                    civic_gene_abbrevations and/or cgi_gene_abbreviations
```
### loh
The *loh* or "loss of heterozygosity" too identifies variants that appear as heterozygotes
in the normal samples, but as homozygotes in the tumor samples. Default settings expect an
allele frequency between 0.3 and 0.7 in the normal samples and exceeded that of 0.8 for
the tumor samples. These values can be adjusted from their defaults with the following
usage options:

```
usage: gemini loh [-h] [--minDP INTEGER] [--minGQ INTEGER] [--maxNorm FLOAT]
                  [--minNorm FLOAT] [--minTumor FLOAT] [--patient STRING]
                  [--samples STRING] [--columns STRING] [--filter STRING]
                  [--purity] [--specific STRING] [--cancers STRING]
                  db

positional arguments:
  db                 The name of the database to be queried.

optional arguments:
  -h, --help         show this help message and exit
  --minDP INTEGER    Minimum depth required in all samples default is 0)
  --minGQ INTEGER    Minimum genotype quality required in all samples (default
                     is 0)
  --maxNorm FLOAT    Specify a maximum normal sample AF to allow (default is
                     0.7)
  --minNorm FLOAT    Specify a minimum normal sample AF to allow (default is
                     0.3)
  --minTumor FLOAT   Specify a minimum AF for tumor samples to require
                     (default is 0.8)
  --patient STRING   Specify a patient to filter (should correspond to a
                     patient_id in ped file)
  --samples STRING   Rather than including all samples, enter a string of
                     comma-separated specified samples to use (default is
                     "All")
  --columns STRING   A comma-separated list of columns that you would like
                     returned (default is "*", which returns every column)
  --filter STRING    Restrictions to apply to variants (SQL syntax)
  --purity           Using purity estimates in cancer manidest, make
                     corrections to AF to be used
  --specific STRING  Search for LOH variants in a single sample compared to
                     the sample(s) that precede it (must specify single sample
                     included among --samples, also --minNorm, --maxNorm will
                     now apply to the preceding sample)
  --cancers STRING   Restrict results to variants/genes associated with
                     specific cancer types by entering a comma-separated
                     string of cancer type abbreviations (see documents for
                     abbreviations) REQUIRES that db include
                     civic_gene_abbrevations and/or cgi_gene_abbreviations
```
### truncal
The *truncal* tool recovers variants that appear to be present in all included tumor
samples, but absent from all normal samples. By default it will require that the
allele frequency of any variant be 0 in the normal samples, but greater than that
in all tumor samples. These requirements can be adjusted with the following usage options:

```
usage: gemini truncal [-h] [--minDP INTEGER] [--minGQ INTEGER]
                      [--maxNorm FLOAT] [--patient STRING] [--samples STRING]
                      [--increase FLOAT] [--columns STRING] [--filter STRING]
                      [--purity] [--somatic_only] [--cancers STRING]
                      db

positional arguments:
  db                The name of the database to be queried.

optional arguments:
  -h, --help        show this help message and exit
  --minDP INTEGER   Minimum depth required in all samples default is 0)
  --minGQ INTEGER   Minimum genotype quality required in all samples (default
                    is 0)
  --maxNorm FLOAT   Optional: specify a maximum normal sample AF to allow
                    (default is 0)
  --patient STRING  Specify a patient to filter (should correspond to a
                    patient_id in ped file)
  --samples STRING  Optional: rather than including all samples, a string of
                    comma-separated specified samples to use (default is
                    "All")
  --increase FLOAT  Optional: add amount to increase truncal AF filter between
                    normal and tumor samples (default is 0)
  --columns STRING  A list of columns that you would like returned (default is
                    "*", which returns every column)
  --filter STRING   Restrictions to apply to variants (SQL syntax)
  --purity          Using purity estimates in ped file, make corrections to AF
                    to be used
  --somatic_only    Only include variants that have been marked as somatic via
                    the set_somatic command
  --cancers STRING  Restrict results to variants/genes associated with
                    specific cancer types by entering a comma-separated string
                    of cancer type abbreviations (see documents for
                    abbreviations) REQUIRES that db include
                    civic_gene_abbrevations and/or cgi_gene_abbreviations
```
### unique
To identify variants that appear to be unique to a sample (or group) or sample(s), the
*unique* tool can be used. By default thsi tool expects the allele frequency of all other
non-specified samples that are included to be 0, while all specified samples have an
allele frequency greater than 0. These parameters can be adjusted with the following
usage options:

```
usage: gemini unique [-h] [--minDP INTEGER] [--minGQ INTEGER]
                     [--specific STRING] [--maxOthers FLOAT]
                     [--patient STRING] [--samples STRING] [--increase FLOAT]
                     [--columns STRING] [--filter STRING] [--purity]
                     [--somatic_only] [--cancers STRING]
                     db

positional arguments:
  db                 The name of the database to be queried.

optional arguments:
  -h, --help         show this help message and exit
  --minDP INTEGER    Minimum depth required in all samples default is 0)
  --minGQ INTEGER    Minimum genotype quality required in all samples (default
                     is 0)
  --specific STRING  Identify unique variants that exist only in samples from
                     this comma-separated list
  --maxOthers FLOAT  Specify a maximum sample AF to allow in other samples
                     (default is 0)
  --patient STRING   Specify a patient to filter (should correspond to a
                     patient_id in ped file)
  --samples STRING   Rather than including all samples in filters, a string of
                     comma-separated specified samples to use (default is
                     "All")
  --increase FLOAT   Add amount to increase AF filter between unique and other
                     samples (default is 0)
  --columns STRING   A list of columns that you would like returned (default
                     is "*", which returns every column)
  --filter STRING    Restrictions to apply to variants (SQL syntax)
  --purity           Using purity estimates in ped file, make corrections to
                     AF to be used
  --somatic_only     Only include variants that have been marked as somatic
                     via the set_somatic command
  --cancers STRING   Restrict results to variants/genes associated with
                     specific cancer types by entering a comma-separated
                     string of cancer type abbreviations (see documents for
                     abbreviations) REQUIRES that db include
                     civic_gene_abbrevations and/or cgi_gene_abbreviations
```
Citation
================
If you use OncoGEMINI in your research, please cite the following manuscript:


Acknowledgements
================
OncoGEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Utah and is led by Tom Nicholas and Aaron Quinlan.  Substantial contributions and discussions 
have been made by Michael Cormier, Brent Pedersen, Yi Qiao, Xiaomeng Huang, and Gabor Marth.

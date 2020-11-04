OncoGEMINI
=============================================================================

Overview
========
OncoGEMINI is an adaptation of [GEMINI](https://github.com/arq5x/gemini) intended for the improved identification of
biologically and clincally relevant tumor variants from multi-sample and longitudinal
tumor sequencing data. Using a GEMINI-compatible database (generated from an annotated 
VCF file), OncoGEMINI is able to filter tumor variants based on included genomic
annotations and various allele frequency signatures. 

![overview](https://github.com/fakedrtom/oncogemini/blob/master/images/overview.png)

Installation
============
To create an `oncogemini` executable, first make sure the proper conda channels are added:
```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```
Then simply install `oncogemini`:
```
conda install -c bioconda oncogemini
```
This will also create executables for `vcfanno` and `vcf2db.py`, which OncoGEMINI is designed to work with.

For all OncoGEMINI scripts and files, clone this repo:
```
git clone https://github.com/fakedrtom/oncogemini.git
```
The `setup.py` can also create an `oncogemini` executable, but will not create executables for `vcfanno`
and `vcf2db.py` like the the conda installer will.
```
python setup.py install
```
Test the executable by running the `master-test.sh` script:
```
cd oncogemini
./master-test.sh
```
This will first create several test OncoGEMINI databases and then run through a series of tests that will
see if basic functionalities of various OncoGEMINI tools and commands are functioning as expected. All tests
that pass will be indicated with an ok or the lack of an error. If all tests pass, the test databases and
temporary files generated throughout the tests will then be removed.

Documentation
================
Since OncoGEMINI retains much of the functionality of GEMINI, it may also be 
helpful to refer to GEMINI's official documentation which can be found [here](http://gemini.readthedocs.org/en/latest/).

VCF Preparation
----------------
Like GEMINI, multi-allelic sites need to be decomposed and normalized using [vt](https://genome.sph.umich.edu/wiki/Vt).
A more thorough explanation and guide for doing this can be found at the [GEMINI documentation](https://gemini.readthedocs.io/en/latest/#new-gemini-workflow).
Provided that vt is available and in your path, the following from the GEMINI docs
should be sufficient for decomposing and normalizing a VCF:
```
# setup
VCF=your.vcf.gz
NORMVCF=your.norm.vcf.gz
REF=your_reference.fasta

# decompose, normalize and annotate VCF with snpEff.
# NOTE: can also swap snpEff with VEP
zless $VCF \
   | sed 's/ID=AD,Number=./ID=AD,Number=R/' \
   | vt decompose -s - \
   | vt normalize -r $REF - \
   | bgzip -c > $NORMVCF
tabix -p vcf $NORMVCF
```
Similarly, it is recommended that a VCF be annotated with either VEP or SnpEff before
additional annotations are included.

OncoGEMINI relies upon VCF annotations for the creation of searchable fields within 
a database. Therefore it is important that a VCF be annotated with all information that
a user desires for filtering. With that in mind, OncoGEMINI was designed to be used 
alongside [vcfanno](https://github.com/brentp/vcfanno) to accomplish all VCF annotation needs. 
Please consult the vcfanno link for details regarding its proper usage, but in short, with
a completed vcfanno configuration file, VCFs can be annotated quite simply:
```
vcfanno vcfanno.config prepared.vcf.gz > annotated.vcf
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
tumor samples with different sampling times), and any sample purity values (optional), if known.
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
vcf2db.py annotated.vcf.gz sample.manifest database.db
```

Usage
----------------
OncoGEMINI utilizes SQL queries in combination with tool commands to search the 
database for variants that match requested filters. Some of the more prominant
tools, and examples for using them, are listed below.

### *query*
For general searches, the *query* tool allows for customization. This is carried over
from GEMINI and further details can be found [here](https://gemini.readthedocs.io/en/latest/content/querying.html#basic-queries).
The *query* command is highly flexible and specific. For example, to search for all
variants on chromosome 13 that have a 'HIGH' impact severity, the following *query*
command would return the chromosome, start and end positions, reference and alternate
alleles and gene for all variants that meet those requirements (if any exist):
```
oncogemini query -q "select chrom, start, end, ref, alt, gene from variants where chrom == 13 and impact_severity == 'HIGH'" database.db
```
### *bottleneck*
The *bottleneck* tool is designed to identify variants whose allele frequencies increase across
sampling timepoints. By default, *bottleneck* will require the slope made by all included allele
frequencies to be greater than 0.05, and the R correlation coefficient for all allele frequencies
to be greater than 0.5. If a normal sample has been included, it will also require that variant
allele frequencies for that sample be 0. Please note that the *bottleneck* tool will also require
that all included tumor samples have a positive (> 0) slope. These and other parameters can be
adjusted with the following usage options:
```
optional arguments:
  --maxNorm FLOAT   Specify a maximum normal sample AF to allow (default is 0)
  --minSlope FLOAT  Minimum slope required for the AFs across samples (default
                    is 0.05)
  --minR FLOAT      Minimum r correlation coefficient required for AFs
                    (default is 0.5)
  --minEnd FLOAT    Minimum AF required of the sample representing the final
                    timepoint (default is 0)
  --endDiff FLOAT   Minimum required AF difference between the samples
                    representing the first and final timepoints (default is 0)
```
For example, to find variants that are increasing in allele frequency across all included
samples, but that exhibit a steeper slope and high correlation coefficient, the following
command could be used:
```
oncogemini bottleneck --minSlope 0.4 --minR 0.8 database.db
```
### *loh*
The *loh* or "loss of heterozygosity" tool identifies variants that appear as heterozygotes
in the normal samples, but as homozygotes in the tumor samples. A normal sample must be
including for the *loh* tool to function properly. Default settings expect an
allele frequency between 0.3 and 0.7 in the normal samples and exceeded that of 0.8 for
the tumor samples. These values can be adjusted from their defaults with the following
usage options:
```
optional arguments:
  --maxNorm FLOAT    Specify a maximum normal sample AF to allow (default is
                     0.7)
  --minNorm FLOAT    Specify a minimum normal sample AF to allow (default is
                     0.3)
  --minTumor FLOAT   Specify a minimum AF for tumor samples to require
                     (default is 0.8)
  --specific STRING  Search for LOH variants in a single sample compared to
                     the sample(s) that precede it (must specify single sample
                     included among --samples, also --minNorm, --maxNorm will
                     now apply to the preceding sample)
```
To more narrowly define heterozygozity in the normal samples and increase the
homozygozity threshold in the tumor samples, the defaults can be changed:
```
oncogemini loh --maxNorm 0.6 --minNorm 0.4 --minTumor 0.9 database.db
```
To identify a loss of heterozygozity variant in a single sample rather than across
all tumor samples compared to the normal samples, the `--specific` parameter can be used.
In this case, the *loh* tool will focus on the specified sample and compare it to the
sample(s) that most immediately precede it, as indicated in the sample manifest's time
column. Heterozygozity in this preceding sample(s) is defined by the `--maxNorm` and
`--minNorm` parameters (or their defaults). For example, if the samples A0, A1, A2, and
A3 (with times indicated as 0, 1, 2, and 3) were loaded into the database, to identify
loss of heterozygozity variants in only A3, the following command is used:
```
oncogemini loh --specific A3 --maxNorm 0.55 --minNorm 0.45 database.db
```
In this example case, the `--maxNorm` and `--minNorm` parameters would be applied to
the A2 sample.

### *truncal*
The *truncal* tool recovers variants that appear to be present in all included tumor
samples, but absent from all normal samples. A normal sample muct be included for
the *truncal* tool to work. By default it will require that the allele frequency of
any variant be 0 in the normal samples, but greater than that in all tumor samples.
These requirements can be adjusted with the following usage options:
```
optional arguments:
  --maxNorm FLOAT   Optional: specify a maximum normal sample AF to allow
                    (default is 0)
  --increase FLOAT  Optional: add amount to increase truncal AF filter between
                    normal and tumor samples (default is 0)
```
Here is a command that would allow for variants with non-zero allele frequencies in
the normal sample(s) and require that the tumor samples have allele frequencies that
are at least 0.2 greater than the maximum allowed allele frequencies in the normal
sample(s):
```
oncogemini truncal --maxNorm 0.05 --increase 0.2 database.db
```
### *unique*
To identify variants that appear to be unique to a sample (or group of samples), the
*unique* tool can be used. By default this tool expects the allele frequency of all other
non-specified samples that are included to be 0, while all specified samples have an
allele frequency greater than 0. These parameters can be adjusted with the following
usage options:
```
optional arguments:
  --specific STRING  Identify unique variants that exist only in samples from
                     this comma-separated list
  --maxOthers FLOAT  Specify a maximum sample AF to allow in other samples
                     (default is 0)
  --increase FLOAT   Add amount to increase AF filter between unique and other
                     samples (default is 0)
```
If the database contains samples B0, B1, and B2, the *unique* tool can identify variants
that are only found in B2:
```
oncogemini unique --specific B2 database.db
```
### Common Parameters
The *bottleneck*, *loh*, *truncal*, and *unique* tools share the following parameters:
```
optional arguments:
  -h, --help         show this help message and exit
  --minDP INTEGER    Minimum depth required in all samples default is 0)
  --minGQ INTEGER    Minimum genotype quality required in all samples (default
                     is 0)
  --patient STRING   Specify a patient to filter (should correspond to a
                     patient_id in ped file)
  --samples STRING   Rather than including all samples, enter a string of
                     comma-separated specified samples to use (default is
                     "All")
  --columns STRING   A comma-separated list of columns that you would like
                     returned (default is "*", which returns every column)
  --filter STRING    Restrictions to apply to variants (SQL syntax)
  --purity           Using purity estimates in cancer manifest file, make
                     corrections to AF to be used
  --somatic_only    Only include variants that have been marked as somatic via
                    the set_somatic command
  --cancers STRING  Restrict results to variants/genes associated with
                    specific cancer types by entering a comma-separated string
                    of cancer type abbreviations (see documents for 
                    abbreviations) REQUIRES that db include
                    civic_gene_abbrevations and/or cgi_gene_abbreviations
```
Of particular note are the `--columns` and `--filter` parameters. With `--columns` the desired
output is specified while `--filter` allows for the listing of variant requirements. For example,
`--columns "chrom, start, end, ref, alt, gene"` and `--filter "impact_severity != 'LOW' and gene ==
'BRCA2'"` will return the chromosome, start and end positions, reference and alternate allele, 
and gene name for any variants that have an impact severity of 'MED' or 'HIGH' and are located 
within the *BRCA2* gene. These are both highly customizable. If `--columns` is not invoked, all 
information for a given variant that is stored in the database will be returned and if `--filter` 
is not used, the variants will not be filtered with any criteria other than those that are built 
into provided tools.

The `--cancers` parameter allows filtered results to be limited to variants in genes with
reported associations with specific cancer types. Currently this is intended to be used alongside
annotations from CIViC and CGI and is not available for use without these annotations (please refer
to the [CRAB](https://github.com/fakedrtom/crab) to include these annotations). For a list of cancer
types and their accepted abbreviations, please refer to [this](https://github.com/fakedrtom/crab/blob/master/cancer_names_abbreviations.txt).

### Somatic Mutations
OncoGEMINI will evaluate all variants within the database and select those that meet specified tool
and annotation filter requirements. Thus, if the VCF used to create the database contained both
germline and somatic mutations, both mutation types would be considered by OncoGEMINI commands. To
focus solely on somatic mutations, it is recommended that the VCF used for the creation of a OncoGEMINI
database be pre-filtered to only include somatic mutations or that somatic mutations be clearly labeled
in the VCF so they are incorporated as a filterable annotation within the database. If that is not
possible, the *set_somatic* tool may be employed which allows for variants within a OncoGEMINI database
to be “flagged” as somatic based on user defined criteria regarding normal and tumor genotypes or sample sequencing
depths and allele frequencies. OncoGEMINI tools may then take advantage of the `--somatic-only`
parameter to restrict variant evaluations to only those variants that have been marked as somatic in
the database by the *set_somatic* tool.

### *set_somatic*
By default the *set_somatic* tool flags variants as somatic if all normal samples provided are genotyped
as homozygous for the reference allele and at least one of the included tumor samples is genotyped as
heterozygous or homozygous for the alternate allele. Users may override these defaults by providing more
detailed requirements regarding allele frequencies, sequencing depths, and alternate read counts in both
the normal and tumor samples, thus allowing more specific designation of variants that should be flagged
as somatic or not. The following parameters are available to *set_somatic* for defining potential somatic
mutations:
```
optional arguments:
  -h, --help            show this help message and exit
  --minDP MINDP         Minimum depth required in all samples (default is 0)
  --minGQ MINGQ         Minimum genotype quality required in all samples
                        (default is 0)
  --normAF NORMAF       The maximum frequency of the alternate allele in the
                        normal sample (default 0).
  --normCount NORMCOUNT
                        The maximum count of the alternate allele in the
                        normal sample (default 0).
  --normDP NORMDP       The minimum depth allowed in the normal sample to
                        believe somatic (default 0).
  --tumAF TUMAF         The minimum frequency of the alternate allele in the
                        tumor sample (default 0).
  --tumCount TUMCOUNT   The minimum count of the alternate allele in the tumor
                        sample (default 0).
  --tumDP TUMDP         The minimum depth allowed in the tumor sample to
                        believe somatic (default 0).
  --dry-run             Don't set the is_somatic flag, just report what
                        _would_ be set. For testing parameters.
```
If none of the additional normal sample parameters are invoked (`--normAF`, `--normCount`, or `--normDP`)
then the default of all normal samples must be genotyped as homozygous reference for the given variant will
be used. Similarly, if none of the additional tumor sample parameters are invoked (`--tumAF`, `--tumCount`,
or `--tumDP`) then the default of at least one tumor sample being genotyped as heterozygous or homozygous for
the alternate allele is used. For example, the following command will use genotype defaults for both the normal
and tumor samples included in the database and only variants that are entirely genotyped as homozygous for the
reference allele in the normal samples and at least one of the included tumor samples is genotyped as heterozygous
or homozygous for the alternate allele will be flagged as somatic by the *set_somatic* tool:

```
oncogemini set_somatic database.db
```
By invoking any of the normal or tumor sample parameters, the genotype defaults can be replaced with more
specific criteria. For example, to require that somatic variants include at least a single tumor sample with
a higher alternate allele frequency (AF >= 0.2), but otherwise keep the genotype defaults for the normal samples,
we can include the `--tumAF` parameter:

```
oncogemini set_somatic --tumAF 0.2 database.db
```
Similarly we can replace the genotype defaults for the normal and tumor samples all at once. With the following
command we can allow that variants be flagged as somatic while exhibiting greater than 0 alternate allele
frequencies in normal samples, while also requiring a specific read depth for all normal samples and create a
minimum AF to be found in at least a single tumor sample. Using *set_somatic* with the following options will mark
variants in the database as somatic if these criteria are met:

```
oncogemini set_somatic --normAF 0.05 --normDP 20 --tumAF 0.2 database.db
```
It is important to note that *set_somatic* is **NOT** a somatic variant caller. However, in the absence of proper somatic
variant calls, the *set_somatic* tool enables users to define criteria that is acceptable to them as being consistent
with their expectations for a somatic variant.

Citation
================
If you use OncoGEMINI in your research, please cite [this manuscript](https://www.biorxiv.org/content/10.1101/2020.03.10.979591v1).


Acknowledgements
================
OncoGEMINI is being developed in the [Quinlan lab](http://quinlanlab.org/) at the University
of Utah and is led by Tom Nicholas.

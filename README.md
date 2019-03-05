Cancer-GEMINI
=============================================================================

Overview
========
Cancer-GEMINI is an adaptation of GEMINI intended for the improved identification of
biologically and clincally relevant tumor variants from multi-sample and longitudinal
tumor sequencing data. Using a GEMINI-compatible database (generated from an annotated 
VCF file), Cancer-GEMINI is able to filter tumor variants based on included genomic
annotations and various allele frequency signatures. 

Cancer-GEMINI filters variants within a GEMINI-compatible database which 
integrates genetic variation alongside genomic annotations 

The intent of ``GEMINI`` (GEnome MINIing) is to provide a simple, flexible, and 
powerful framework for exploring genetic variation for personal and medical genetics.
GEMINI is unique in that it integrates genetic variation (from VCF files) with
a wealth of genome annotations into a unified database framework. Using this
integrated database as the analysis framework, we aim to leverage the expressive 
power of SQL for data analysis, while attempting to overcome the fundamental 
challenges associated with using databases for very large
(e.g. 1,000,000 variants times 1,000 samples yields one billion genotypes)
datasets. In addition, by defining sample relationships with a PED file, GEMINI allows
one to explore and test for variants that meet specific inheritance models (e.g., 
recessive, dominant, etc.).


Documentation
================
The official documentation is here:

Since Cancer-GEMINI retains much of the functionality of GEMINI, it may also be 
helpful to refer to GEMINI's official documentation which can be found [here](http://gemini.readthedocs.org/en/latest/)


Citation
================
If you use Cancer-GEMINI in your research, please cite the following manuscript:


Acknowledgements
================
Cancer-GEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Utah and is led by Tom Nicholas and Aaron Quinlan.  Substantial contributions and discussions 
have been made by Brent Pedersen, Yi Qiao, Xiaomeng Huang, and Gabor Marth.



Installation
============
Install ``GEMINI`` using the automated installation script, `gemini_install.py`. This
script installs GEMINI along with required python libraries, third party tools and data 
files used for variant annotation. The installation documentation contains additional 
details on installed files and tools.

http://gemini.readthedocs.org/en/latest/content/installation.html

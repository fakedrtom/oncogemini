Cancer-GEMINI
=============================================================================

Overview
========
Cancer-GEMINI is an adaptation of [GEMINI](https://github.com/arq5x/gemini) intended for the improved identification of
biologically and clincally relevant tumor variants from multi-sample and longitudinal
tumor sequencing data. Using a GEMINI-compatible database (generated from an annotated 
VCF file), Cancer-GEMINI is able to filter tumor variants based on included genomic
annotations and various allele frequency signatures. 


Installation
============


Documentation
================
The official documentation is here.

Since Cancer-GEMINI retains much of the functionality of GEMINI, it may also be 
helpful to refer to GEMINI's official documentation which can be found [here](http://gemini.readthedocs.org/en/latest/).

VCF Preparation
----------------
Like GEMINI, multi-allelic sites need to be decomposed and normalized using [vt](https://genome.sph.umich.edu/wiki/Vt).
A more thorough explanation and guide for doing this can be found at the [GEMINI documentation](https://gemini.readthedocs.io/en/latest/#new-gemini-workflow).

Cancer-GEMINI relies upon VCF annotations for the creation of searchable fields within 
a database. Therefore it is important that a VCF be annotated with all information that
a user desires filtering. Cancer-GEMINI was designed to be used alongside [vcfanno](https://github.com/brentp/vcfanno) to 
accomplish all VCF annotation needs. 
 
Citation
================
If you use Cancer-GEMINI in your research, please cite the following manuscript:


Acknowledgements
================
Cancer-GEMINI is being developed in the Quinlan lab (quinlanlab.org) at the University
of Utah and is led by Tom Nicholas and Aaron Quinlan.  Substantial contributions and discussions 
have been made by Brent Pedersen, Yi Qiao, Xiaomeng Huang, and Gabor Marth.

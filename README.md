Copy Number Assessment of Population Scale Array-CGH (canny)
======================================================

Description
-----------

This repository contains software that was used to analyze 2135 Agilent 1M CGH arrays designed to assess copy number variation genotypes in Phase 3 samples from the 1000 Genomes Project.

Required third-party resources 
------------------------------

A number of third party software packages and modules are required by these programs:

*           R:  http://www.r-project.org/, with LPCM and ggplot2 packages
*    samtools:  https://github.com/samtools/htslib (also contains tabix, bgzip) 
*    vcftools:  http://vcftools.sourceforge.net/
*        perl:  Statistics::R package, http://search.cpan.org/~fangly/Statistics-R-0.32/

In addition, you will need:

* tabix indexed file with log2 ratios for each probe and sample
  -obtainable here ()
* sample list in same order as probe file
  -obtainable here ()

Example workflow
----------------
An example workflow would be as follows:

Input is VCF:

~~~
canny.pl \
--input_filename=test.vcf.gz \
--ratio_filename=allData.filtered.gc.txt.corr3.allCols.dat.gz \ 
--sample_filename=samples.txt \
--plot_results \
--bandwidth=0.05 \
--threshold=0.3 \
--output_filename=test_arraygenotypes.vcf.gz`
~~~

Input is Region:

~~~
canny.pl \
--input_region=7:143911018-144034817 \
--ratio_filename=allData.filtered.gc.txt.corr3.allCols.dat.gz \
--sample_filename=samples.txt \
--plot_results \
--bandwidth=0.05 \
--threshold=0.3 \
--output_filename=test_arraygenotypes.vcf.gz
~~~

Contact
-------
Questions: Please contact Ryan Mills at remills@umich.edu
05/05/2014

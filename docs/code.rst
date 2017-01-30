Appendix 1 - SDGS Validation
*******************************

Comparison with the Genome in a Bottle Reference Sample
=======================================================

This script carries out a comparison of sequencing results to the reference calls for NA12878 (the flagship genome for the Genome in a Bottle Consortium).
Sequencing of a reference sample is now required to meet best practice guidelines and most NHS laboratories are currently using this genome.

In order for the variants to be filtered to the regions of interest an intersection of the BED files for the panel being validated and the GIAB truth regions must be generated. Any regions outside the truth regions cannot be confirmed as homozygous reference and these regions cannot be annotated as true or false calls. 

Generate BED Intersects
-----------------------
.. automethod:: validation.giab_comparison.generate_bed_intersects

Generate Whole BED
+++++++++++++++++++
.. automethod:: validation.giab_comparison.generate_whole_bed

Generate Remainder
+++++++++++++++++++
.. automethod:: validation.giab_comparison.generate_remainder

Compare VCFs
-----------------------
The BED files generated are used to filter the VCF files for comparison using bcftools isec.

Prepare VCF
+++++++++++++++++++
.. automethod:: validation.giab_comparison.prepare_vcf

Intersect VCFs
+++++++++++++++++++
.. automethod:: validation.giab_comparison.bcftools_isec

Annotate results
+++++++++++++++++++
The files produced contain the variants that are shared between the two VCF files and those that are unique to each data set (i.e. false positives and negatives). The false calls are used to generate the sensitivity and MCC values for each panel.
False calls are annotated with the coverage at that position to determine whether there is evidence of the missed variant for false negatives and if there is a coverage related reason for false positives or mismatching genotypes.

False Positives
###################
.. automethod:: validation.giab_comparison.annotate_false_pos

False Negatives
###################
.. automethod:: validation.giab_comparison.annotate_false_negs

Mismatching Genotypes
######################
.. automethod:: validation.giab_comparison.check_genotype

Check Coverage
###################
.. automethod:: validation.giab_comparison.get_coverage

Calculating Remainder
######################

.. automethod:: validation.giab_comparison.remainder_size


Generation of Coverage Graphs for Validation Samples Sequenced
===============================================================

This script creates graphs showing a breakdown of coverage across the validation patients (-d) for the regions in the BED files listed (-b).

Two sets of graphs are generated: proportion of regions covered at pre-defined depths and a breakdown of each exon in each gene to help highlight any areas that may cause problems when the panel goes into service. The graphs are saved in the specified folder (-o) and named either by the panel or gene.

Generating Dictionaries
-------------------------
The regions of interest and corresponding coverage values are used to populate a series of dictionaries. The information within the dictionaries can then be accessed quickly when generating the graphs.

Generate BED Dictionary
++++++++++++++++++++++++

.. automethod:: validation.coverage_plots.generate_bed_dict

Generate Gene Dictionary
++++++++++++++++++++++++++

.. automethod:: validation.coverage_plots.generate_gene_dict

Getting Coverage Values
++++++++++++++++++++++++++

.. automethod:: validation.coverage_plots.get_coverage


Generating Plots
-------------------------

Sub Panel Graphs
++++++++++++++++++

.. automethod:: validation.coverage_plots.plot_region_coverage

Gene Graphs
+++++++++++++++++

.. automethod:: validation.coverage_plots.generate_gene_plots


Calculating the Number of Unique Variants Analysed During the Validation
===========================================================================

This script calculates the number of unique variants identified that are shared, unique to the MiSeq data and unique to the HiSeq data. The output is a series of JSON files detailing the variants and counts.

.. automethod:: validation.UniqueVariants.variant_counts

.. automethod:: validation.UniqueVariants.generate_v_dict

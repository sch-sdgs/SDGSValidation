SDGS Validation Usage
*********************

This page describes how to run each of the relevant scripts on the command line. The arguments and their format etc.
is given along with an example command.

giab_comparison.py
==================

Purpose
-------

This script is designed to compare a reference genome "truth set" to the results of NGS sequencing in house. The script
creates a BED file of regions that overlap between the panel of interest and the truth regions for the reference genome.
It then compares the VCF produced by the pipeline to the truth VCF for the reference genome, using the BED file
generated previously.

The variants are then summarised and a results file is produced. This lists the number of matching variants as well as
the details for any variants that are false positive, negative or mismatching genotype. The script also uses the BAM
file from the pipeline to give coverage levels of the false calls to help determine the reason for the discrepancy.

Arguments
---------

**-o**

Full path to the output folder for the results and any interim files generated. If the folder does not already exist,
the script will create the directory.

**-p**

The sample name as it appears in the VCF. This is used to identify the sample specific variant information, for example
the genotype, from the VCF. Therefore, it is important that this parameter matches exactly to the VCF file column
header.

**-b**

The prefix associated with any BED files that are to be used, for example, NGD\_. This is used to find all BED files
that should be used in the evaluation of the sequencing results. This argument should also contain the path to the BED
files.

If you wanted to include all the current NGD BED files used in the pipeline this argument would be:

.. code-block:: bash

    -b /results/Analysis/MiSeq/MasterBED/NGD_

**-bam**

The full path to the BAM files for the relevant pipeline run. This is the BAM file used to generate the VCF that will
be compared. The BAM file is used to generate coverage information for the variants that are not concordant.

**-v**

The VCF file to be compared to the truth set.

**-rv**

The path to the reference genome VCF. This is the VCF of the "truth variants" that should be compared to during the
analysis. *The default for this argument is
'/results/Analysis/HiSeq_validation/giab/giab-NA12878/truth_small_variants.decomposed.normalised.vcf.gz' which is the
path for the GIAB truth set.* Therefore, if this is the reference sample used, the argument can be omitted.

**-rs**

The identifier used in the reference VCF. This is similar to the *-p* argument but for the "truth set" VCF. *The default
for this argument is "INTEGRATION" which is the identifier used in the GIAB sample truth set.* Therefore, if this
reference sample has been used, the argument can be omitted.

Execution
---------

To run the code, you will also need the following libraries:

* argparse
* pybedtools
* pysam
* json2html

You will also need v2.17.0 of BEDtools and v? sambamba. The path to the bin directory should be in the PATH variable:

.. code-block:: bash

    export PATH=${PATH}:/results/Pipeline/program/bedtools-2.17.0/bin

An example command is given below for comparison of the GIAB sample, anything within the <> symbols should be replaced
with the true values.

.. code-block:: bash

    python giab_comparison.py -o *<output_dir>* -p *<sample_id>* -b *<bed_file_prefix>* -bam *<path_to_BAM_file>*
    -v *<path_to_VCF>*

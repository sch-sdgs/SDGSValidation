import argparse
import glob
import os
import subprocess
from pybedtools import BedTool
from pysam import VariantFile
from math import sqrt
import json
from json2html import *

"""
.. module:: giab_comparison
   :platform: Unix
   :synopsis: Comparison of sequencing results to the reference calls for NA12878 (the flagship genome for the Genome
   in a Bottle Consortium).

.. moduleauthor:: Natalie Groves


"""

def generate_whole_bed(truth_regions, bedtool_list, out_dir):
    """
    Generate broad BED file for the regions that overlap with the truth set

    :param truth_regions: GIAB truth regions
    :type truth_regions: BedTool
    :param bedtool_list: List of BED files in the panel (as BEDTools)
    :type bedtool_list: List of BedTool objects
    :param out_dir: location for whole bed file
    :type out_dir: String
    :return: File path to the final BED file
    :rtype: String
    """
    whole_region = truth_regions.intersect(bedtool_list)
    whole_region_sorted = whole_region.sort()
    whole_region_merged = whole_region_sorted.merge()
    print('out_dir = ' + out_dir)
    whole_bed = out_dir + '/whole.bed'
    whole_region_merged.moveto(whole_bed)
    return whole_bed

def generate_remainder(whole_bed, out_dir, bed_list):
    """
    Calculate the remaining regions that are not included in the truth set

    :param whole_bed: Path to the truth regions for the whole panel
    :type whole_bed: String
    :param out_dir: Prefix used for all the bed files
    :type out_dir: String
    :param bed_list: List of all the bed files for that panel
    :type bed_list: List of String
    :return: BEDTool containing any regions that are completely missing from the truth regions
    :rtype: BedTool
    """
    try:
        whole_truth = BedTool(whole_bed)
        whole_truth.saveas()
        whole = BedTool()

        for bed in bed_list:
            print(bed)
            tool = BedTool(bed)
            tool.saveas()
            if bed == bed_list[0]:
                whole = tool
            else:
                whole = whole.cat(tool)
                whole.saveas()

        whole_sorted = whole.sort()
        whole_merged = whole_sorted.merge()
        whole_merged.saveas()

        remainder = whole_merged.subtract(whole_truth)
        remainder.moveto(out_dir + '/remainder.bed')
        missing_regions = whole_merged.subtract(whole_truth, A=True)
    except UnicodeDecodeError:
        missing_regions = None
    return missing_regions

def generate_bed_intersects(bed_prefix, directory):
    """
    Creates the intersected BED file for the broad panel and each sub panel associated with the given abbreviation

    :param bed_prefix: The prefix given to each of the BED files within the panel
    :type bed_prefix: String
    :param directory: Location of pipeline output
    :type directory: String
    :return: Dictionary containing the abbreviations for each of the one-based bed files
    :rtype: Dictionary
    """
    print('Getting BED files.')
    path = bed_prefix + "*.bed"
    bed_files = glob.glob(path)
    print(bed_files)

    if len(bed_files) == 0:
        print('No BED files found with that prefix')
        exit(1)

    bed_dict = {}
    bedtool_list = []
    truth_regions = BedTool('/results/Analysis/projects/HeredCancer/CysticValidation/HuRef/truth_regions.bed')

    print('Generating truth regions.')
    for f in bed_files:
        name = os.path.basename(f)
        no_header = directory + "/" + name.replace('.bed', '_noheader.bed')
        one_based = directory + "/" + name.replace('.bed', '_truth_regions_1based.bed')
        truth_regions_panel = directory + "/" + name.replace('.bed', '_truth_regions.bed')

        #Create BED file without header for intersect
        command = "grep -i -v start " + f + " > " + no_header
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print('Error executing command: ' + str(e.returncode))
            exit(1)

        #Intersect BEDs
        panel_bed = BedTool(no_header)

        intersect = truth_regions.intersect(panel_bed)
        intersect.saveas(truth_regions_panel)
        bedtool_list.append(intersect)

        #Create one-based BED file for use with bcftools
        command = "awk '{print($1\"\t\"$2+1\"\t\"$3)}' " + truth_regions_panel + " > " + one_based
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print('Error executing command: ' + str(e.returncode))
            exit(1)

        #find the bed abbreviation from the file on the server to name the result files
        command_grep = "grep " + name + " /results/Analysis/MiSeq/MasterBED/abbreviated_bed_names.txt | cut -f2"
        command_length = "awk '{SUM += $3-$2} END {print SUM}' " + truth_regions_panel
        try:
            abv = subprocess.check_output(command_grep, shell=True)
            bed_length = subprocess.check_output(command_length, shell=True)
            try:
                int(bed_length)
            except ValueError:
                bed_length = 0
            bed_dict[one_based] = {'abv':abv.strip(), 'length':bed_length}
        except subprocess.CalledProcessError as e:
            print('Error executing command: ' + str(e.returncode))
            exit(1)
        os.remove(no_header)

    whole_bed = generate_whole_bed(truth_regions, bedtool_list, directory)
    missing_regions = generate_remainder(whole_bed, directory, bed_files)
    if missing_regions:
        missing_regions.moveto(directory + '/missing_regions.bed')

    print('BED files produced correctly.')
    return bed_dict, whole_bed

def prepare_vcf(vcf):
    """
    Vcf must be decomposed, normalised and zipped and indexed before it can be used with bcftools

    :param vcf: The file path to the original vcf
    :type vcf: String
    :return: File path to the decomposed and zipped vcf
    :rtype: String
    """
    print('Preparing vcf.')
    decomposed = vcf.replace('.vcf', '.decomposed.vcf')
    normalised = vcf.replace('.vcf', '.decomposed.normalised.vcf')
    normalised_zipped = normalised + '.gz'

    try:
        command = '/results/Pipeline/program/vt/vt decompose ' +  vcf + ' -o ' + decomposed
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)

    try:
        command = '/results/Pipeline/program/vt/vt normalize ' + decomposed + ' -r ' + \
                  '/results/Pipeline/program/GATK_resource_bundle/ucsc.hg19.nohap.masked.fasta -o ' + normalised
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)

    try:
        command = 'bgzip ' + normalised
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    try:
        command = 'tabix ' + normalised_zipped
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    print('vcf decomposed and zipped successfully.')
    return normalised_zipped

def get_coverage(whole_bed, directory, sample, bam):
    """
    Coverage at all positions is calculated.

    This is then used for coverage analysis and to determine read depth at any false negative sites

    :param whole_bed: BED file generated from all BED files in panel
    :type whole_bed: String
    :param directory: Location of patient results
    :type directory: String
    :param sample: sample ID
    :type sample: String
    :param bam: Path to BAM file to be used in coverage calculation
    :type bam: String
    :return: Filename for coverage stats
    :rtype: String
    """
    print('Generating coverage stats.')
    if directory.endswith('.txt'):
        out = directory
        directory = os.path.dirname(directory)
    else:
        out = directory + '/whole_bed_coverage.txt'
    print(whole_bed)
    command = '/results/Pipeline/program/sambamba/build/sambamba depth base --min-coverage=0 -q29 -m -L ' + whole_bed + \
              ' ' + bam + ' > ' + out + '.tmp'
    print(command)
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)
    print('Sambamba complete.')
    #issue with sambamba that leaves out regions that have 0 coverage - intersect regions to find missing and add
    # them to the file at coverage 0
    temp_bed = out.replace('.txt', '.bed.tmp')
    command = 'awk \'{print($1"\\t"$2"\\t"$2+1"\\t"$3)}\' ' + out + '.tmp | grep -v "COV" > ' + temp_bed
    print(command)
    try:
        subprocess.check_call(command, shell=True)
        print('BED coordinates extracted.')
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)

    coverage_bed = BedTool(temp_bed)
    print('BED tool created')
    whole_bedtool = BedTool(whole_bed)
    print('Intersecting')
    missing_regions = whole_bedtool.intersect(coverage_bed, v=True)
    missing_file = directory + 'regions_missing'
    missing_regions.moveto(missing_file)
    print('Generating file')
    command = '''while read i; do start=`echo "$i"|cut -f2`; end=`echo "$i"|cut -f3`; chr=`echo "$i"|cut -f1`; end_true=`echo "${end} - 1" | bc`; for j in $(seq $start $end_true); do new_end=`echo -e "${j} + 1" | bc`; echo -e "$chr\\t${j}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t''' + sample + '''";done;done < ''' + missing_file + '> ' + directory + '/to_add'
    print(command)
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)

    command = 'cat ' + out + '.tmp ' + directory + '/to_add > ' + out
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        exit(1)
    print('fix complete.')
    return out

def annotate_false_negs(folder, ref_sample, coverage_file):
    """
    Get information for any false negative results.

    Returns basic variant info plus quality, genotype, coverage (total, ref base and alt base if appropriate)

    :param folder: Folder containing output from bcftools isec
    :type folder: String
    :param ref_sample: Sample number for reference vcf
    :type ref_sample: String
    :param coverage_file: File containing per base coverage for the truth_regions panel
    :type coverage_file: String
    :return: List of variant dictionaries containing information on false negatives
    :rtype: List
    """
    false_negs = VariantFile(folder + '/0000.vcf')
    num_neg = len(list(false_negs.fetch()))
    print(num_neg)

    variants = {'indels':[],'no_coverage':[],'evidence_of_alt':[],'false_neg':[]}
    count=0
    if num_neg > 0:
        print('false negatives')
        for rec in false_negs.fetch():
            chrom = rec.contig
            pos = int(rec.pos)
            ref = rec.alleles[0]
            alt = rec.alleles[1]
            qual = rec.qual
            genotype = rec.samples[ref_sample]['GT']
            count+=1
            if len(rec.alleles[0]) == 1 and len(rec.alleles[1]) == 1:
                search = '\'' + rec.contig + '\s' + str(rec.pos - 1) + '\''
                command = 'grep ' + search + ' ' + coverage_file
                try:
                    line = subprocess.check_output(command, shell=True)
                except subprocess.CalledProcessError as e:
                    print(command)
                    print('Error executing command: ' + str(e.returncode))
                    exit(1)

                if line == '':
                    variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                               'GT':genotype, 'coverage':{'total':'no coverage information', 'ref':'N/A', 'alt':'N/A'}}
                    no_cov = variants['no_coverage']
                    no_cov.append(variant)
                    variants['no_coverage'] = no_cov
                else:
                    line.strip('\n')
                    bases = {'A': 3, 'C': 4, 'G': 5, 'T': 6}
                    fields = line.split()
                    cov = fields[2]
                    ref_cov = fields[bases[rec.alleles[0]]]
                    alt_cov = fields[bases[rec.alleles[1]]]
                    variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                               'GT':genotype, 'coverage':{'total':cov, 'ref':ref_cov, 'alt':alt_cov}}

                    if int(cov) == 0:
                        no_cov = variants['no_coverage']
                        no_cov.append(variant)
                        variants['no_coverage'] = no_cov
                    elif int(alt_cov) != 0:
                        ev_alt = variants['evidence_of_alt']
                        ev_alt.append(variant)
                        variants['evidence_of_alt'] = ev_alt
                    else:
                        fn = variants['false_neg']
                        fn.append(variant)
                        variants['false_neg'] = fn
            else:
                variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                            'GT':genotype, 'coverage':{'total':'indel: no coverage could be obtained', 'ref':'N/A',
                                                        'alt':'N/A'}}
                indels = variants['indels']
                indels.append(variant)
                variants['indels'] = indels

    else:
        print('no false negatives')
    print("false_negatives=" + str(count))
    return variants

def annotate_false_pos(folder, coverage_file, sample):
    """
    Get information for any false positive results

    Returns basic variant info plus quality, genotype, coverage (total, ref base and alt base if appropriate)

    :param folder: Folder containing output from bcftools isec
    :type folder: String
    :param coverage_file: File containing per base coverage for the truth_regions panel
    :type coverage_file: String
    :param sample: Container ID used in vcf file
    :type sample: String
    :return: List of variant dictionaries containing information on false negatives
    :rtype: List
    """
    false_pos = VariantFile(folder + '/0001.vcf')
    num_pos = len(list(false_pos.fetch()))
    print(num_pos)

    variants = []

    if num_pos > 0:
        print('false positives')
        for rec in false_pos.fetch():
            chrom = rec.contig
            pos = int(rec.pos)
            ref = rec.alleles[0]
            alt = rec.alleles[1]
            qual = rec.qual
            genotype = rec.samples[sample]['GT']
            if 'AD' in rec.samples[sample].keys():
                allelic_depth = rec.samples[sample]['AD']
            elif 'NV' in rec.samples[sample].keys():
                allelic_depth = rec.samples[sample]['NV']
            else:
                allelic_depth = 'N/A'
            if 'DP' in rec.samples[sample].keys():
                total_depth = rec.samples[sample]['DP']
            elif 'NR' in rec.samples[sample].keys():
                total_depth = rec.samples[sample]['NR']
            else:
                total_depth = 0
            if len(rec.alleles[0]) == 1 and len(rec.alleles[1]) == 1:
                search = '\'' + rec.contig + '\s' + str(rec.pos - 1) + '\''
                command = 'grep ' + search + ' ' + coverage_file
                try:
                    line = subprocess.check_output(command, shell=True)
                except subprocess.CalledProcessError as e:
                    print('Error executing command: ' + str(e.returncode))
                    exit(1)
                if line == '':
                    variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                               'GT':genotype, 'vcf_depth':{'DP':total_depth, 'AD':allelic_depth},
                               'coverage':{'total':'no coverage information', 'ref':'N/A', 'alt':'N/A'}}
                else:
                    bases = {'A': 3, 'C': 4, 'G': 5, 'T': 6}
                    fields = line.split()
                    cov = fields[2]
                    ref_cov = fields[bases[rec.alleles[0]]]
                    alt_cov = fields[bases[rec.alleles[1]]]
                    variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                               'GT':genotype, 'vcf_depth':{'DP':total_depth, 'AD':allelic_depth},
                               'coverage':{'total':cov, 'ref':ref_cov, 'alt':alt_cov}}
            else:
                variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                            'GT':genotype, 'vcf_depth':{'DP':total_depth, 'AD':allelic_depth},
                           'coverage':{'total':'indel: no coverage could be obtained', 'ref':'N/A', 'alt':'N/A'}}
            variants.append(variant)
    else:
        print('no false positives')

    return variants

def check_genotype(folder, sample, ref_sample, coverage_file):
    """
    Compares the genotype for all shared variants

    The number of matching variants are counted and those that do not match are annotated with basic variant info plus
    quality, genotype, coverage (total, ref base and alt base if appropriate)

    :param folder: Location of results from the NGS analysis pipeline
    :type folder: String
    :param sample: Sample number (used in vcf file)
    :type sample: String
    :param ref_sample: Sample number for reference vcf
    :type ref_sample: String
    :param coverage_file: File containing coverage information for each position in the panel
    :type coverage_file: String
    :return: Number of matching variants
    :rtype: Int
    :return: List of variant dictionaries with detailed information for mismatching genotypes
    :rtype: List
    """
    shared_giab = VariantFile(folder + '/0002.vcf')
    shared_patient = VariantFile(folder + '/0003.vcf')

    variants = []

    vars_giab = {}
    for rec in shared_giab.fetch():
        chrom = rec.contig
        pos = rec.pos
        alleles = rec.alleles
        if chrom not in vars_giab:
            vars_giab[chrom] = {}
        if pos not in vars_giab[chrom]:
            vars_giab[chrom][pos] = {}
        if alleles not in vars_giab[chrom][pos]:
            vars_giab[chrom][pos][alleles] = rec.samples[ref_sample]['GT']

    matching = 0
    for rec in shared_patient.fetch():
        chrom = rec.contig
        pos = rec.pos
        alleles = rec.alleles
        if 'AD' in rec.samples[sample].keys():
            allelic_depth = rec.samples[sample]['AD']
        else:
            allelic_depth = 'N/A'
        if 'DP' in rec.samples[sample].keys():
            total_depth = rec.samples[sample]['DP']
        elif 'NR' in rec.samples[sample].keys():
            total_depth = rec.samples[sample]['NR']
        else:
            total_depth = 0
        giab_genotype = vars_giab[chrom][pos][alleles]
        if rec.samples[sample]['GT'] == giab_genotype:
            matching += 1
        elif (rec.samples[sample]['GT'][0] is None or rec.samples[sample]['GT'][0] == 1) and rec.samples[sample]['GT'][
            0] == giab_genotype[1] and rec.samples[sample]['GT'][1] == giab_genotype[0]:
            matching += 1
        elif rec.samples[sample]['GT'][0] == 0 and rec.samples[sample]['GT'][1] == 1 and giab_genotype[0] == 1 and giab_genotype[1] == 0:
            matching += 1
        elif rec.samples[sample]['GT'][0] == 1 and rec.samples[sample]['GT'][1] == 0 and giab_genotype[0] == 0 and giab_genotype[1] == 1:
            matching += 1
        else:
            if len(rec.alleles[0]) == 1 and len(rec.alleles[1]) == 1:
                search = '\'' + rec.contig + '\s' + str(rec.pos - 1) + '\''
                command = 'grep ' + search + ' ' + coverage_file
                try:
                    line = subprocess.check_output(command, shell=True)
                except subprocess.CalledProcessError as e:
                    print('Error executing command: ' + str(e.returncode))
                    exit(1)
                if line == '':
                    variant = {'chrom': chrom, 'pos': pos, 'ref': alleles[0], 'alt': alleles[1], 'QUAL': rec.qual,
                               'GT': {"sample": rec.samples[sample]['GT'], 'GIAB': giab_genotype},
                               'vcf_depth': {'DP': total_depth, 'AD': allelic_depth},
                               'coverage':{'total':'no coverage information', 'ref':'N/A', 'alt':'N/A'}}
                else:
                    bases = {'A': 3, 'C': 4, 'G': 5, 'T': 6}
                    fields = line.split()
                    cov = fields[2]
                    ref_cov = fields[bases[rec.alleles[0]]]
                    alt_cov = fields[bases[rec.alleles[1]]]
                    variant = {'chrom':chrom, 'pos':pos, 'ref':alleles[0], 'alt':alleles[1], 'QUAL':rec.qual,
                                    'GT':{"sample":rec.samples[sample]['GT'], 'GIAB':giab_genotype},
                                    'vcf_depth':{'DP':total_depth, 'AD':allelic_depth},
                                    'coverage':{'total':cov, 'ref':ref_cov, 'alt':alt_cov}}
            else:
                variant = {'chrom': chrom, 'pos': pos, 'ref': alleles[0], 'alt': alleles[1], 'QUAL': rec.qual,
                           'GT': {"sample": rec.samples[sample]['GT'], 'GIAB': giab_genotype},
                           'vcf_depth': {'DP': total_depth, 'AD': allelic_depth},
                           'coverage': {'total': 'indel: no coverage could be obtained', 'ref': 'N/A', 'alt': 'N/A'}}
            variants.append(variant)
    print(str(matching) + ' matching variants')

    return matching, variants

def remainder_size(bed_file):
    """
    Calculate the number of bases in the small panel not included in the truth regions.

    The amount of the panel that overlaps with the truth regions is calculated. This gives an idea of how well the
    panel is represented in the reference sample and whether it is an accurate reflection of the accuracy of the
    process.

    :param bed_file: File path to the specific panel bed file
    :type bed_file: String
    :return: Total length of regions in remainder BED file
    :rtype: Int
    """
    print('Calculating remainder')
    print(bed_file)
    original_bed = '/results/Analysis/projects/GIAB/merged_fastqs/' + os.path.basename(bed_file.replace('_truth_regions_1based', ''))
    print(original_bed)
    truth_region_bed = bed_file.replace('_1based', '')
    print(truth_region_bed)

    truth_region_tool = BedTool(truth_region_bed)

    original_tool = BedTool(original_bed)

    remainder_name = bed_file.replace('_truth_regions_1based', '_remainder')
    remainder_bed = original_tool.subtract(truth_region_tool)
    remainder_bed.saveas(remainder_name)

    f = open(remainder_name, 'r')
    regions = [line.strip('\n') for line in f.readlines()]
    f.close()

    total_length = 0

    for region in regions:
        fields = region.split('\t')
        start = int(fields[1])
        end = int(fields[2])

        length = end - start
        total_length += length

    return total_length

def bcftools_isec(file_prefix, decomposed_zipped, bed_dict, bam, reference_vcf, reference_sample, whole_bed, out_dir):
    """
    Intersect the two vcfs and limit to the truth regions and panel BED file.

    The method counts the number of false positives and false negatives and checks the genotype of all of the matching
    variants.

    Any false negatives are investigated in terms of depth

    :param file_prefix: Prefix given to all files during the NGS pipeline (i.e. worklist-patient)
    :type file_prefix: String
    :param decomposed_zipped: File path for the decomposed and zipped vcf
    :type decomposed_zipped: String
    :param bed_dict: Dictionary containing the abbreviations for each of the BED files - to be used as folder names
    :type bed_dict: Dictionary
    :param bam: Path to BAM file to be used in coverage calculation
    :type bam: String
    :param reference_vcf: Path to the reference genome vcf
    :type reference_vcf: String
    :param reference_sample: Identifier in reference sample vcf
    :type reference_sample: String
    :param whole_bed: Path to combined BED file
    :type whole_bed: String
    :param out_dir: Directory for results files
    :type out_dir: String
    :return: Analysis of variant comparison
    :rtype: Dictionary
    """
    print('Comparing vcfs.')
    sample_split = file_prefix.split('-')
    sample = sample_split[1] + '-' + sample_split[2]

    results = out_dir + '/giab_results'
    try:
        os.mkdir(results)
    except OSError as e:
        print(results)
        print(e.errno)
        if e.errno != 17:
            exit(1)

    coverage_file = get_coverage(whole_bed, results, sample, bam)

    path = out_dir + '/*truth_regions_1based.bed'
    print(path)
    bed_files = glob.glob(path)

    all_results = {}
    print(bed_files)
    for f in bed_files:
        abv = bed_dict[f]['abv']
        print(abv)
        folder = results + '/' + abv
        try:
            os.mkdir(folder)
        except OSError as e:
            if e.errno != 17:
                print(folder)
                print(e.errno)
                exit(1)

        if os.path.getsize(f) > 0:

            command = '/results/Pipeline/program/bcftools-1.3.1/bcftools isec -R ' + f + ' -p ' + folder + \
                      ' ' + reference_vcf +' ' + decomposed_zipped
            try:
                subprocess.check_call(command, shell=True)
            except subprocess.CalledProcessError as e:
                print('Error executing command: ' + str(e.returncode))
                exit(1)

            print('Checking for false calls.')
            false_negs_ann = annotate_false_negs(folder, reference_sample, coverage_file)
            false_pos_ann = annotate_false_pos(folder, coverage_file, sample)
            print('Checking genotype for shared calls.')
            num_matching, genotypes = check_genotype(folder, sample, reference_sample, coverage_file)

            false_negs_indels = len(false_negs_ann['indels'])
            false_negs_cov = len(false_negs_ann['no_coverage'])
            false_negs_ev_alt = len(false_negs_ann['evidence_of_alt'])
            false_neg_oth = len(false_negs_ann['false_neg'])
            false_negs = false_neg_oth + false_negs_cov + false_negs_ev_alt + false_negs_indels
            print("false neg total=" + str(false_negs))
            false_pos = len(false_pos_ann)
            num_mismatch = len(genotypes)
            true_positives = num_matching + num_mismatch

            total_bases = int(bed_dict[f]['length'])
            true_negatives = total_bases - (false_negs + false_pos + num_mismatch + num_matching)

            if true_negatives == 0:
                print('ERROR: Coverage file empty')
                exit(1)

            try:
                sensitivity = (num_matching + num_mismatch) / float((true_positives + false_pos)) * 100
            except ZeroDivisionError:
                sensitivity = 0.0

            total_pos = true_positives + false_pos
            actual_pos = true_positives + false_negs
            total_negs = true_negatives + false_negs
            actual_negs = true_negatives + false_pos

            if total_pos == 0 or actual_pos == 0 or total_negs == 0 or actual_negs == 0:
                denominator = 1
            else:
                denominator = total_pos * total_negs * actual_pos * actual_negs

            mcc = (true_positives * true_negatives - false_pos * false_negs) / \
                  sqrt(denominator)

            remainder_length = remainder_size(f)
            percent_covered = float(total_bases) / (total_bases + remainder_length) * 100

            out = {'false_negative': false_negs_ann, 'false_positive': false_pos_ann,
                   'mismatching_genotype': genotypes, 'matching_variants': num_matching,
                   'num_true_negatives': true_negatives, 'sensitivity': sensitivity, 'MCC': mcc,
                   'small_panel_remainder_length': remainder_length, 'percent_small_panel_covered': percent_covered,
                   'num_false_positive': false_pos, 'num_false_negative':{'indel':false_negs_indels,
                                                                          'no_coverage':false_negs_cov,
                                                                          'ev_of_alt':false_negs_ev_alt,
                                                                          'false_neg':false_neg_oth,
                                                                          'total':false_negs},
                   'num_mismatching_genotype': num_mismatch}
        else:
            total_bases = int(bed_dict[f]['length'])
            print(total_bases)

            original_bed = '/results/Analysis/projects/GIAB/merged_fastqs/' + os.path.basename(f).replace('_truth_regions_1based', '')
            print(original_bed)

            f = open(original_bed, 'r')
            regions = [line.strip('\n') for line in f.readlines()]
            f.close()

            remainder_length = 0

            for region in regions:
                if region.startswith('#'):
                    pass
                else:
                    fields = region.split('\t')
                    start = int(fields[1])
                    end = int(fields[2])

                    length = end - start
                    remainder_length += length

            print(remainder_length)
            percent_covered = float(total_bases) / (total_bases + remainder_length) * 100
            print(percent_covered)

            out = {'false_negative': {}, 'false_positive': {},
                   'mismatching_genotype': {}, 'matching_variants': 0,
                   'num_true_negatives': 0, 'sensitivity': 0, 'MCC': 0,
                   'small_panel_remainder_length': remainder_length, 'percent_small_panel_covered': percent_covered,
                   'num_false_positive':0, 'num_false_negative':{'indel':0, 'no_coverage':0, 'ev_of_alt':0, 'false_neg':0, 'total':0},
                   'num_mismatching_genotype':0}
        all_results[abv] = out

    return all_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='Location of output', required=True)
    parser.add_argument('-p', help='Prefix given to all results files (i.e. worklist-patient)', required=True)
    parser.add_argument('-b', help='Prefix associated with panel BED files (i.e. NGD_)', required=True)
    parser.add_argument('-bam', help='BAM file for run', required=True)
    parser.add_argument('-v', help='VCF to be compared', required=True)
    parser.add_argument('-rv', help='Path to reference genome vcf', default='/results/Analysis/HiSeq_validation/giab/giab-NA12878/truth_small_variants.decomposed.normalised.vcf.gz')
    parser.add_argument('-rs', help='Identifier used in reference vcf', default='INTEGRATION')

    args = parser.parse_args()

    directory=args.o
    if directory.endswith('/'):
        directory = os.path.dirname(directory)
    vcf = args.v
    file_prefix=args.p
    bed_prefix=args.b
    bam = args.bam
    reference_vcf = args.rv
    reference_sample = args.rs

    command = 'export PATH=${PATH}:/results/Pipeline/program/bedtools-2.17.0/bin'
    print(command)
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(command)
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    bed_dict, whole_bed = generate_bed_intersects(bed_prefix, directory)

    decomposed_zipped = prepare_vcf(vcf)

    results = bcftools_isec(file_prefix, decomposed_zipped, bed_dict, bam, reference_vcf, reference_sample, whole_bed, directory)

    f = open(directory+'/giab_summary.txt', 'w')
    j = json.dumps(results, indent=4)
    print >> f, j
    f.close()

    f = open(directory + '/giab_summary.html', 'w')
    f.write(json2html.convert(json=j))
    f.close()





if __name__ == '__main__':
    main()
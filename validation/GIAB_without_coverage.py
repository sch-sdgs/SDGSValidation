import argparse
import glob
import os
import subprocess
from pybedtools import BedTool
from pysam import VariantFile
from math import sqrt
import json

def generate_whole_bed(truth_regions, bedtool_list, bed_prefix):
    """
    Generate broad BED file for the regions that overlap with the truth set
    :param truth_regions: GIAB truth regions
    :param bedtool_list: list of BED files in the panel (as BEDTools)
    :param bed_prefix: Prefix used for all BED files in the panel
    :return: filepath to the final BED file
    """
    whole_region = truth_regions.intersect(bedtool_list)
    whole_region_sorted = whole_region.sort()
    whole_region_merged = whole_region_sorted.merge()
    whole_bed = '/results/Analysis/MiSeq/MasterBED/GIAB/' + bed_prefix + '.whole.bed'
    whole_region_merged.moveto(whole_bed)
    return whole_bed

def generate_remainder(whole_bed, bed_prefix, bed_list):
    """
    Calculate the remaining regions that are not included in the truth set
    :param whole_bed: path to the truth regions for the whole panel
    :param bed_prefix: prefix used for all the bed files
    :param bed_list: list of all the bed files for that panel
    :return: BEDTool containing any regions that are completely missing from the truth regions
    """

    whole_truth = BedTool(whole_bed)
    whole_truth.saveas()
    whole = BedTool()

    for bed in bed_list:
        print bed
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
    remainder.moveto('/results/Analysis/MiSeq/MasterBED/GIAB/' + bed_prefix + '.remainder.bed')
    missing_regions = whole_merged.subtract(whole_truth, A=True)
    return missing_regions

def generate_bed_intersects(bed_prefix, directory):
    """
    Creates the intersected BED file for the broad panel and each sub panel associated with the given abbreviation
    :param bed_prefix: The prefix given to each of the BED files within the panel
    :param directory: location of pipeline output
    :return: dictionary containing the abbreviations for each of the one-based bed files
    """
    print 'Getting BED files.'
    path = '/results/Analysis/MiSeq/MasterBED/' + bed_prefix + "*"
    bed_files = glob.glob(path)

    if len(bed_files) == 0:
        print 'No BED files found with that prefix'
        exit(1)

    bed_dict = {}
    bedtool_list = []
    truth_regions = BedTool('/results/Analysis/MiSeq/MasterBED/GIAB/truth_regions.bed')

    print 'Generating truth regions.'
    for f in bed_files:
        name = os.path.basename(f)
        no_header = '/results/Analysis/MiSeq/MasterBED/GIAB/' + name.replace('.bed', '_noheader.bed')
        one_based = '/results/Analysis/MiSeq/MasterBED/GIAB/' + name.replace('.bed', '_truth_regions_1based.bed')
        truth_regions_panel = '/results/Analysis/MiSeq/MasterBED/GIAB/' + name.replace('.bed', '_truth_regions.bed')

        #Create BED file without header for intersect
        command = "grep -i -v start " + f + " > " + no_header
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print 'Error executing command: ' + str(e.returncode)
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
            print 'Error executing command: ' + str(e.returncode)
            exit(1)

        #find the bed abbreviation from the file on the server to name the result files
        command_grep = "grep " + name + " /results/Analysis/MiSeq/MasterBED/abbreviated_bed_names.txt | cut -f2"
        command_length = "awk '{SUM += $3-$2} END {print SUM}' " + no_header
        try:
            abv = subprocess.check_output(command_grep, shell=True)
            bed_length = subprocess.check_output(command_length, shell=True)
            bed_dict[one_based] = {'abv':abv.strip(), 'length':bed_length}
        except subprocess.CalledProcessError as e:
            print 'Error executing command: ' + str(e.returncode)
            exit(1)
        os.remove(no_header)

    whole_bed = generate_whole_bed(truth_regions, bedtool_list, bed_prefix)
    missing_regions = generate_remainder(whole_bed, bed_prefix, bed_files)
    missing_regions.moveto(directory + '/missing_regions.bed')

    print 'BED files produced correctly.'
    return bed_dict

def prepare_vcf(vcf):
    """
    vcf must be decomposed, zipped and indexed before it can be used with bcftools
    :param vcf: The vcf file to be compared
    :return: filepath to the decomposed and zipped vcf
    """
    print 'Preparing vcf.'

    decomposed = vcf.replace('.vcf', '.decomposed.vcf')
    normalised = vcf.replace('.vcf', '.decomposed.normailsed.vcf')
    normalised_zipped = normalised + '.gz'
    print vcf
    try:
        command = '/results/Pipeline/program/vt/vt decompose ' +  vcf + ' -o ' + decomposed
        print command
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print 'Error executing command:' + str(e.returncode)
        exit(1)

    try:
        command = '/results/Pipeline/program/vt/vt normalize ' + decomposed + ' -r ' + \
                  '/results/Pipeline/program/GATK_resource_bundle/ucsc.hg19.nohap.masked.fasta -o ' + normalised
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print 'Error executing command:' + str(e.returncode)
        exit(1)

    try:
        command = 'bgzip ' + normalised
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print 'Error executing command: ' + str(e.returncode)
        exit(1)

    try:
        command = 'tabix ' + normalised_zipped
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print 'Error executing command: ' + str(e.returncode)
        exit(1)

    print 'vcf decomposed and zipped successfully.'
    return normalised_zipped

def annotate_false_negs(folder):
    """
    Get information for any false negative results - returns basic variant info plus quality, genotype, coverage
    (total, ref base and alt base if appropriate)
    :param folder: Folder containing output from bcftools isec
    :return: array of variant dictionaries containing information on false negatives
    """
    false_negs = VariantFile(folder + '/0000.vcf')
    num_neg = len(list(false_negs.fetch()))
    print num_neg

    variants = []

    if num_neg > 0:
        print 'false negatives'
        for rec in false_negs.fetch():
            chrom = rec.contig
            pos = int(rec.pos)
            ref = rec.alleles[0]
            alt = rec.alleles[1]
            qual = rec.qual
            genotype = rec.samples['INTEGRATION']['GT']

            variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                       'GT':genotype}

            variants.append(variant)
    else:
        print 'no false negatives'

    return variants

def annotate_false_pos(folder, sample):
    """
    Get information for any false positive results - returns basic variant info plus quality, genotype, coverage
    (total, ref base and alt base if appropriate)
    :param folder: Folder containing output from bcftools isec
    :param sample: container ID used in vcf file
    :return: array of variant dictionaries containing information on false negatives
    """
    false_pos = VariantFile(folder + '/0001.vcf')
    num_pos = len(list(false_pos.fetch()))
    print num_pos

    variants = []

    if num_pos > 0:
        print 'false positives'
        for rec in false_pos.fetch():
            chrom = rec.contig
            pos = int(rec.pos)
            ref = rec.alleles[0]
            alt = rec.alleles[1]
            qual = rec.qual
            genotype = rec.samples[sample]['GT']
            if 'AD' in rec.samples[sample].keys():
                allelic_depth = rec.samples[sample]['AD']
            else:
                allelic_depth = 'N/A'
            total_depth = rec.samples[sample]['DP']
            variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt, 'QUAL':qual,
                       'GT':genotype, 'vcf_depth':{'DP':total_depth, 'AD':allelic_depth}}

            variants.append(variant)
    else:
        print 'no false positives'

    return variants

def check_genotype(folder, sample):
    """
    Compares the genotype for all shared variants
    :param folder: location of results from the NGS analysis pipeline
    :param sample:  sample number (used in vcf file)
    :return: dictionary of number of matching variants and detailed information for any with mismatching genotypes
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
            vars_giab[chrom][pos][alleles] = rec.samples['INTEGRATION']['GT']

    matching = 0
    for rec in shared_patient.fetch():
        chrom = rec.contig
        pos = rec.pos
        alleles = rec.alleles
        if 'AD' in rec.samples[sample].keys():
            allelic_depth = rec.samples[sample]['AD']
        else:
            allelic_depth = 'N/A'
        total_depth = rec.samples[sample]['DP']
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

            variant = {'chrom': chrom, 'pos': pos, 'ref': alleles[0], 'alt': alleles[1], 'QUAL': rec.qual,
                       'GT': {sample: rec.samples[sample]['GT'], 'GIAB': giab_genotype},
                       'vcf_depth': {'DP': total_depth, 'AD': allelic_depth}}

            variants.append(variant)
    print str(matching) + ' matching variants'
    results = {'matching':matching, 'mismatching':variants}
    print results
    return results

def remainder_size(bed_prefix):
    """
    Calculate the number of bases in the broad panel not included in the truth regions
    :param bed_prefix: prefix used for all bed files in the panel (used to create remainder file name
    :return:total length of regions in bed file
    """
    remainder  = '/results/Analysis/MiSeq/MasterBED/GIAB/' + bed_prefix + '.remainder.bed'
    f = open(remainder, 'r')
    regions = [line.strip('\n') for line in f.readlines()]

    total_length = 0

    for region in regions:
        fields = region.split('\t')
        start = int(fields[1])
        end = int(fields[2])

        length = end - start
        total_length += length

    return total_length

def bcftools_isec(sample, decomposed_zipped, bed_prefix, bed_dict):
    """
    Intersect the two vcfs and limit to the truth regions and panel BED file.
    The method counts the number of false positives and false negatives and checks the genotype of all of the matching
    variants.
    Any false negatives are investigated in terms of depth
    :param sample: sample name as appears in the vcf
    :param decomposed_zipped: Filepath for the decomposed and zipped vcf
    :param bed_prefix: Prefix for the BED files in the panel
    :param bed_dict: Dictionary containing the abbreviations for each of the BED files - to be used as folder names
    :return: Analysis of variant comparison
    """
    print 'Comparing vcfs.'

    directory = os.path.dirname(decomposed_zipped)
    results = directory + '/giab_results'
    try:
        os.mkdir(results)
    except OSError as e:
        print results
        print e.errno
        if e.errno != 17:
            exit(1)

    path = '/results/Analysis/MiSeq/MasterBED/GIAB/' + bed_prefix + '*truth_regions_1based.bed'
    bed_files = glob.glob(path)

    all_results = {}

    for f in bed_files:
        abv = bed_dict[f]['abv']
        print abv
        folder = results + '/' + abv
        try:
            os.mkdir(folder)
        except OSError as e:
            if e.errno != 17:
                print folder
                print e.errno
                exit(1)

        command = '/results/Pipeline/program/bcftools-1.3.1/bcftools isec -R ' + f + ' -p ' + folder + \
                  ' /results/Analysis/HiSeq_validation/giab/giab-NA12878/truth_small_variants.decomposed.vcf.gz ' + \
                  decomposed_zipped
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print 'Error executing command: ' + str(e.returncode)
            exit(1)

        print 'Checking for false calls.'
        false_negs_ann = annotate_false_negs(folder)
        false_pos_ann = annotate_false_pos(folder, sample)
        print 'Checking genotype for shared calls.'
        genotypes = check_genotype(folder, sample)

        false_negs = len(false_negs_ann)
        false_pos = len(false_pos_ann)
        num_matching = genotypes['matching']
        num_mismatch = len(genotypes['mismatching'])
        true_positives = num_matching + num_mismatch

        total_bases = int(bed_dict[f]['length'])
        true_negatives = total_bases - (false_negs + false_pos + num_mismatch + num_matching)

        if true_negatives == 0:
            print 'ERROR: Coverage file empty'
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

        remainder_length = remainder_size(bed_prefix)
        percent_covered = float(total_bases) / (total_bases + remainder_length) * 100

        out = {'false_negative': false_negs_ann, 'false_positive': false_pos_ann,
               'mismatching_genotype': genotypes['mismatching'], 'matching_variants': genotypes['matching'],
               'num_true_negatives': true_negatives, 'sensitivity': sensitivity, 'MCC': mcc,
               'broad_panel_remainder_length': remainder_length, 'percent_broad_panel_covered': percent_covered}

        all_results[abv] = out

    return all_results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', help='vcf to be compared', required=True)
    parser.add_argument('-s', help='sample id (as appears in vcf)', required=True)
    parser.add_argument('-b', help='Prefix associated with panel BED files (i.e. NGD_)', required=True)


    args = parser.parse_args()

    vcf=args.v
    directory = os.path.dirname(vcf)
    sample=args.s
    bed_prefix=args.b

    bed_dict = generate_bed_intersects(bed_prefix, directory)

    decomposed_zipped = prepare_vcf(vcf)

    results = bcftools_isec(sample, decomposed_zipped, bed_prefix, bed_dict)

    f = open(directory+'/giab_summary.txt', 'w')
    j = json.dumps(results, indent=4)
    print >> f, j
    f.close()





main()
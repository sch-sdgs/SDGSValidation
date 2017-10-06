from giab_comparison import get_coverage, annotate_false_negs, annotate_false_pos, check_genotype, prepare_vcf
import argparse
import os
import subprocess
import json
from pybedtools import BedTool
from math import sqrt

def isec_vcf(out=None,v1=None,s1=None,v1_prep=None,v2=None,s2=None,bed=None,bam=None,ref=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='Location of output')
    parser.add_argument('-v1', help='VCF 1')
    parser.add_argument('-s1', help='Sample name (as appears in vcf 1)')
    parser.add_argument('-v1_prep', help='Does the reference VCF need to be prepared', default=True)
    parser.add_argument('-v2', help='VCF 2')
    parser.add_argument('-s2', help='Sample name (as appears in vcf 2), if omitted it will be the same as s1', default=None)
    parser.add_argument('-b', help='BED file to filter to')
    parser.add_argument('-bam', help='BAM file for vcf 2 sample to generate coverage stats from')
    parser.add_argument('-ref', help="Path to the Reference Genome used", default='/results/Pipeline/program/GATK_resource_bundle/ucsc.hg19.nohap.masked.fasta')

    args = parser.parse_args()

    if not v1:
        v1 = args.v1
    if not s1:
        s1 = args.s1
    if not v2:
        v2 = args.v2
    if not s2 and not args.s2:
        s2 = s1
    elif not s2:
        s2 = args.s2
    if not bed:
        bed = args.b
    if not bam:
        bam = args.bam
    if not v1_prep:
        v1_prep = args.v1_prep

    if not out:
        out = args.o

    if not ref:
        ref = args.ref

    if out.endswith('/'):
        out_dir = os.path.dirname(out)
    else:
        out_dir = out

    if not os.path.exists(out_dir):
        print("Output directory does not exist")
        exit(1)

    no_header = os.path.basename(bed).replace('.bed', '_noheader.bed')
    no_header = out_dir + '/' + no_header
    command = "grep -i -v start " + bed + " > " + no_header
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    if os.path.basename(no_header) == "IEM_all_panels_header_noheader.bed":
        #IEM name column causes bedtools error
        tool = BedTool(no_header)
        sorted_tool = tool.sort()
        merged_tool = sorted_tool.merge()
        merged_tool.saveas(no_header.replace('.bed', '.merged.bed'))
        no_header = no_header.replace('.bed', '.merged.bed')

    print(no_header)

    #prepare vcfs
    if v1_prep != "False":
        v1_decomposed = prepare_vcf(v1, ref)
    else:
        print('no prep')
        v1_decomposed = v1
    if os.path.exists(v2.replace('.vcf', '.decomposed.normalised.vcf.gz')):
        v2_decomposed = v2.replace('.vcf', '.decomposed.normalised.vcf.gz')
    else:
        v2_decomposed = prepare_vcf(v2, ref)

    #isec
    command = '/results/Pipeline/program/bcftools-1.3.1/bcftools isec -R ' + \
                bed + ' -p ' + \
                out_dir + ' ' + \
                v1_decomposed + ' ' + \
                v2_decomposed
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    #annotate
    coverage = get_coverage(no_header, out_dir, s2, bam)
    false_negs_ann = annotate_false_negs(out_dir, s1, coverage)
    false_pos_ann = annotate_false_pos(out_dir, coverage, s2)
    matching, mismatching_genotype = check_genotype(out_dir, s2, s1, coverage)
    false_negs_indels = len(false_negs_ann['indels'])
    false_negs_cov = len(false_negs_ann['no_coverage'])
    false_negs_ev_alt = len(false_negs_ann['evidence_of_alt'])
    false_neg_oth = len(false_negs_ann['false_neg'])
    false_negs = false_neg_oth + false_negs_cov + false_negs_cov + false_negs_ev_alt + false_negs_indels
    false_pos = len(false_pos_ann)
    num_mismatch = len(mismatching_genotype)
    true_positives = matching + num_mismatch

    try:
        sensitivity = matching / float((true_positives + false_pos)) * 100
    except ZeroDivisionError:
        sensitivity = 0.0

    f = open(coverage, 'r')
    lines = f.readlines()
    total_bases = len(lines) -1
    f.close()
    true_negatives = total_bases - (false_negs + false_pos + num_mismatch + matching)

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

    out = {'result':{'false_negative': false_negs_ann, 'false_positive': false_pos_ann,
           'mismatching_genotype': mismatching_genotype, 'matching_variants': matching,
           'num_true_negatives': true_negatives, 'sensitivity': sensitivity, 'MCC': mcc,
           'small_panel_remainder_length': "Not calculated", 'percent_small_panel_covered': "Not calculated",
           'num_false_positive': false_pos, 'num_false_negative': {'indel': false_negs_indels,
                                                                   'no_coverage': false_negs_cov,
                                                                   'ev_of_alt': false_negs_ev_alt,
                                                                   'false_neg': false_neg_oth,
                                                                   'total': false_negs},
           'num_mismatching_genotype': num_mismatch}}

    f = open(out_dir + '/summary.json', 'w')
    j = json.dumps(out, indent=4)
    print >> f, j
    f.close()

if __name__ == '__main__':
    isec_vcf()
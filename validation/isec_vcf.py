from giab_comparison import get_coverage, annotate_false_negs, annotate_false_pos, check_genotype, prepare_vcf
import argparse
import os
import subprocess
import json
from pybedtools import BedTool

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='Location of output', required=True)
    parser.add_argument('-v1', help='VCF 1', required=True)
    parser.add_argument('-s1', help='Sample name (as appears in vcf 1)', required=True)
    parser.add_argument('-v2', help='VCF 2', required=True)
    parser.add_argument('-s2', help='Sample name (as appears in vcf 2), if omitted it will be the same as s1', default=None)
    parser.add_argument('-b', help='BED file to filter to', required=True)
    parser.add_argument('-bam', help='BAM file for vcf 2 sample to generate coverage stats from', required=True)

    args = parser.parse_args()

    v1 = args.v1
    s1 = args.s1
    v2 = args.v2
    if not args.s2:
        s2 = s1
    else:
        s2 = args.s2
    bed = args.b
    bam = args.bam

    if args.o.endswith('/'):
        out_dir = os.path.dirname(args.o)
    else:
        out_dir = args.o

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
    if os.path.exists(v1.replace('.vcf', '.decomposed.normalised.vcf.gz')):
        v1_decomposed = v1.replace('.vcf', '.decomposed.normalised.vcf.gz')
    else:
        v1_decomposed = prepare_vcf(v1)
    if os.path.exists(v2.replace('.vcf', '.decomposed.normalised.vcf.gz')):
        v2_decomposed = v2.replace('.vcf', '.decomposed.normalised.vcf.gz')
    else:
        v2_decomposed = prepare_vcf(v2)

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
    false_negs = annotate_false_negs(out_dir, s1, coverage)
    false_pos = annotate_false_pos(out_dir, coverage, s2)
    matching, mismatching_genotype = check_genotype(out_dir, s2, s1, coverage)

    out = {'false_negative': false_negs, 'num_false_neg':len(false_negs),
           'false_positive': false_pos, 'num_false_pos':len(false_pos),
           'mismatching_genotype': mismatching_genotype, 'num_mismatching':len(mismatching_genotype),
           'matching_variants': matching}

    f = open(out_dir + '/summary.json', 'w')
    j = json.dumps(out, indent=4)
    print >> f, j
    f.close()

if __name__ == '__main__':
    main()
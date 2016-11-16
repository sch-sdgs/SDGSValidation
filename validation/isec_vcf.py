from giab_comparison import annotate_false_negs, annotate_false_pos, check_genotype, prepare_vcf
import argparse
import os
import subprocess
import json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='Location of output', required=True)
    parser.add_argument('-v1', help='VCF 1', required=True)
    parser.add_argument('-s1', help='Sample name (as appears in vcf 1)', required=True)
    parser.add_argument('-v2', help='VCF 2', required=True)
    parser.add_argument('-s2', help='Sample name (as appears in vcf 2), if omitted it will be the same as s1', default=None)
    parser.add_argument('-b', help='BED file to filter to', required=True)
    parser.add_argument('-c', help='Coverage file (_coverage_depth_bases_full_small_panel.bed)', required=True)

    args = parser.parse_args()

    v1 = args.v1
    s1 = args.s1
    v2 = args.v2
    if not args.s2:
        s2 = s1
    else:
        s2 = args.s2
    bed = args.b
    coverage = args.c

    if args.o.endswith('/'):
        out_dir = os.path.dirname(args.o)
    else:
        out_dir = args.o

    if not os.path.exists(out_dir):
        print("Output directory does not exist")
        exit(1)

    #prepare vcfs
    v1_decomposed = prepare_vcf(v1)
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
    false_negs = annotate_false_negs(out_dir, s1, coverage)
    false_pos = annotate_false_pos(out_dir, coverage, s2)
    matching, mismatching_genotype = check_genotype(out_dir, s2, s1, coverage)

    out = {'false_negative': false_negs, 'false_positive': false_pos,
           'mismatching_genotype': mismatching_genotype, 'matching_variants': matching}

    f = open(out_dir + '/summary.json', 'w')
    j = json.dumps(out, indent=4)
    print >> f, j
    f.close()

if __name__ == '__main__':
    main()
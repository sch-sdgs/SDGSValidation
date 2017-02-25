import argparse
from pysam import VariantFile
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='Results VCF to be compared', required=True)
    parser.add_argument('-bed', help='The reference BED file', required=True)
    parser.add_argument('-s', help='Sample ID in VCF', required=True)
    parser.add_argument('-out', help='The folder to putt results files', required=True)

    args = parser.parse_args()

    if args.out.endswith('/'):
        out_dir = args.out
    else:
        out_dir = args.out + '/'

    sample = args.s

    vcf_file = args.vcf
    bed = args.bed

    f = open(bed, 'r')
    regions = [line.strip('\n') for line in f.readlines()]
    f.close()

    variants = {}
    for region in regions:
        if region.startswith('#'):
            continue
        chrom, start, end, name = region.split('\t')
        pos, ref, alt = name.split(':')
        if chrom not in variants:
            variants[chrom] = {pos:{(ref, alt):False,}}
        elif pos not in variants[chrom]:
            variants[chrom][pos] = {(ref, alt):False,}
        else:
            variants[chrom][pos][(ref, alt)] = False

    vcf = VariantFile(vcf_file)
    false_pos = []
    false_neg = []
    true_pos = []
    for v in vcf.fetch():
        chrom = v.contig
        pos = str(v.pos)
        ref = v.alleles[0]
        alt = v.alleles[1]
        qual = v.qual
        genotype = v.samples[sample]['GT']
        if 'AD' in v.samples[sample].keys():
            allelic_depth = v.samples[sample]['AD']
        elif 'NV' in v.samples[sample].keys():
            allelic_depth = v.samples[sample]['NV']
        else:
            allelic_depth = 'N/A'
        if 'DP' in v.samples[sample].keys():
            total_depth = v.samples[sample]['DP']
        elif 'NR' in v.samples[sample].keys():
            total_depth = v.samples[sample]['NR']
        else:
            total_depth = 0
        if pos in variants[chrom].keys():
            if (ref,alt) in variants[chrom][pos].keys():
                variant = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'QUAL': qual,
                           'GT': genotype, 'vcf_depth': {'DP': total_depth, 'AD': allelic_depth},
                           'coverage': {'total': 'no coverage information', 'ref': 'N/A', 'alt': 'N/A'}}
                true_pos.append(variant)
                variants[chrom][pos][(ref, alt)] = True
            else:
                variant = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'QUAL': qual,
                           'GT': genotype, 'vcf_depth': {'DP': total_depth, 'AD': allelic_depth},
                           'coverage': {'total': 'no coverage information', 'ref': 'N/A', 'alt': 'N/A'}}
                false_pos.append(variant)
        else:
            variant = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'QUAL': qual,
                       'GT': genotype, 'vcf_depth': {'DP': total_depth, 'AD': allelic_depth},
                       'coverage': {'total': 'no coverage information', 'ref': 'N/A', 'alt': 'N/A'}}
            false_pos.append(variant)

    for chrom in variants.keys():
        for pos in variants[chrom].keys():
            for v in variants[chrom][pos].keys():
                if not variants[chrom][pos][v]:
                    variant = {'chrom': chrom, 'pos': pos, 'ref': v[0], 'alt': v[1], 'QUAL': 0,
                               'GT': (0,0),
                               'coverage': {'total': 'no coverage information', 'ref': 'N/A', 'alt': 'N/A'}}
                    false_neg.append(variant)

    out = {'false_negative': {'indels':[],'no_coverage':[],'evidence_of_alt':[],'false_neg':false_neg}, 'false_positive': false_pos,
           'mismatching_genotype': [], 'matching_variants': len(true_pos),
           'num_true_negatives': 0, 'sensitivity': 0, 'MCC': 0,
           'small_panel_remainder_length': 0, 'percent_small_panel_covered': 0,
           'num_false_positive': len(false_pos), 'num_false_negative': {'indel': 0,
                                                                   'no_coverage': 0,
                                                                   'ev_of_alt': 0,
                                                                   'false_neg': 0,
                                                                   'total': len(false_neg)},
           'num_mismatching_genotype': 0}

    all_results = {sample:out}
    f = open(out_dir + sample + '_summary.json', 'w')
    j = json.dumps(all_results, indent=4)
    print >> f, j
    f.close()

if __name__ == "__main__":
    main()
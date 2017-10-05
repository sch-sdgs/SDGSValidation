import argparse
import json

def txt_file():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--filelist")

    args = arg_parser.parse_args()

    f = open(args.filelist, 'r')
    files = [line.strip('\n') for line in f.readlines()]

    unique_variants = []

    for v_file in files:
        f = open(v_file, 'r')
        variants = [line.strip('\n') for line in f.readlines()]
        for var in variants:
            col = var.split('\t')
            chrom = col[0]
            pos = col[1]
            ref = col[2]
            alt = col[3]

            variant = {'chrom':chrom, 'pos':pos, 'ref':ref, 'alt':alt}
            if variant not in unique_variants:
                unique_variants.append(variant)


    print(len(unique_variants))

def v_list():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--file")
    arg_parser.add_argument("--out")

    args = arg_parser.parse_args()

    f = open(args.file, 'r')
    variants = [line.strip('\n') for line in f.readlines()]

    unique_variants = []
    multiple_alleles = 0
    v_dict = {}

    for v in variants:
        fields = v.split('\t')
        chrom = fields[0]
        pos = fields[1]
        ref = fields[2]
        alt = fields[3]

        if chrom not in v_dict:
            v_dict[chrom] = {}
            v_dict[chrom][pos] = {}
            v_dict[chrom][pos][ref] = alt

            unique_variants.append(v)
        elif pos not in v_dict[chrom]:
            v_dict[chrom][pos] = {}
            v_dict[chrom][pos][ref] = alt
            unique_variants.append(v)
        elif v_dict[chrom][pos][ref] != alt:
            multiple_alleles += 1
            unique_variants.append(v)

    f = open(args.out, 'w')
    f.write('\n'.join(unique_variants))
    f.close()
    print(len(unique_variants))
    print(multiple_alleles)

def generate_v_dict(variants):
    """
    Generate a Dictionary of variants from a text file.

    Text files from the pipeline can be in one of two formats. The method determines the format from the header and
    populates the dictionary accordingly. The dictionary keys are a string of chr:pos:ref and the value is a list of
    all the alternate alleles in the file.

    :param variants: List of lines from the variant file
    :type variants: List
    :return: Dictionary of positions with the alternate alleles seen in the file
    :rtype: Dictionary
    """
    v_dict = {}
    old = False
    vcf = False
    for line in variants:
        if line.startswith('#CHROM'):
            vcf = True
            continue
        elif line.startswith('#'):
            continue
        elif line.startswith('Variant'):
            old = True
            continue
        elif line.startswith('chrom'):
            old = False
            continue
        elif vcf:
            v_fields = line.split('\t')
            chrom = v_fields[0]
            pos = v_fields[1]
            ref = v_fields[3]
            alt = v_fields[4]
        elif old:
            v_fields = line.split('\t')
            chrom = v_fields[2]
            pos = v_fields[3]
            ref = v_fields[5]
            alt = v_fields[6]
        else:
            v_fields = line.split('\t')
            chrom = v_fields[0]
            pos = v_fields[1]
            ref = v_fields[2]
            alt = v_fields[3]


        key = chrom + ':' + pos + ':' + ref
        if key in v_dict:
            alts = v_dict[key]
            alts.append(alt)
            v_dict[key] = alts
        else:
            v_dict[key] = [alt,]

    return v_dict

def dict_from_v_list(var_list):
    """
    Creates dictionaries for each sample that ahs had Sanger confirms detailing concordant variants for use in count

    :param var_list: filename for list of concordant variants
    :return: dictionary of patients and corresponding variants
    """
    results = {}
    f = open(var_list)
    lines = [line.strip('\n').strip('\r') for line in f.readlines()]
    f.close()

    for line in lines:
        if line.startswith("S number"):
            continue
        fields = line.split('\t')
        sample = fields[0]
        print(sample)
        if sample in results.keys():
            v_dict = results[sample]
        else:
            v_dict = {}
        chrom = fields[1]
        pos = fields[2]
        ref = fields[3]
        alt = fields[4]

        key = chrom + ':' + pos + ':' + ref
        if key in v_dict:
            alts = v_dict[key]
            alts.append(alt)
            v_dict[key] = alts
        else:
            v_dict[key] = [alt, ]

        results[sample] = v_dict
    return results

def variant_counts():
    """
    Compare two variant files for a number of patients to calculate the number of unique variants and amount of overlap.

    The script calculates the number of shared variants, unique to HiSeq and unique to MiSeq for each patient.
    A count of the overall number of unique variants is also calculated as this must reach at least 60 for the
    validation of a new panel.

    :param --files: Path to a file containing the paths to the variant files. Must be formatted Sample,MiSeq_file,HiSeq_file using commas as separators
    :type --files: String
    :param --out: The path to the folder for the output of the script
    :type --out: String
    :return: A number of json files are saved in the specified output directory: One containing the counts for each patient and the total, and one for each set of variants. The variant jsons include the basic information for each variant (chrom, pos, ref, alt).
    :rtype: JSON files
    """
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("--files")
    arg_parser.add_argument("--out")

    args = arg_parser.parse_args()

    total_unique_shared = {}
    count_unique_shared = 0
    total_unique_miseq = {}
    count_unique_miseq = 0
    total_unique_hiseq = {}
    count_unique_hiseq = 0
    sample_breakdown = {}

    f = open(args.files, 'r')
    samples = [line.strip('\n') for line in f.readlines()]
    f.close()

    for sample in samples:
        fields = sample.split(',')
        s = fields[0]
        print(s)
        if s == "sanger_variants":
            results = dict_from_v_list(fields[1])
            print(results.keys())
            for v in results.keys():
                miseq_var = results[v]
                hiseq_var = results[v]

                sample_shared = 0
                sample_miseq = 0
                sample_hiseq = 0

                for var in miseq_var.keys():
                    for alt in miseq_var[var]:
                        try:
                            if alt in hiseq_var[var]:
                                sample_shared += 1
                                if var not in total_unique_shared:
                                    total_unique_shared[var] = {}
                                    total_unique_shared[var][alt] = 1
                                    count_unique_shared += 1
                                elif alt not in total_unique_shared[var]:
                                    total_unique_shared[var][alt] = 1
                                    count_unique_shared += 1
                                else:
                                    count = total_unique_shared[var][alt]
                                    count += 1
                                    total_unique_shared[var][alt] = count
                            else:
                                sample_miseq += 1
                                if var not in total_unique_miseq:
                                    total_unique_miseq[var] = {}
                                    total_unique_miseq[var][alt] = 1
                                    count_unique_miseq += 1
                                elif alt not in total_unique_miseq[var]:
                                    total_unique_miseq[var][alt] = 1
                                    count_unique_miseq += 1
                                else:
                                    count = total_unique_miseq[var][alt]
                                    count += 1
                                    total_unique_miseq[var][alt] = count

                        except KeyError:
                            sample_miseq += 1
                            if var not in total_unique_miseq:
                                total_unique_miseq[var] = {}
                                total_unique_miseq[var][alt] = 1
                                count_unique_miseq += 1
                            elif alt not in total_unique_miseq[var]:
                                total_unique_miseq[var][alt] = 1
                                count_unique_miseq += 1
                            else:
                                count = total_unique_miseq[var][alt]
                                count += 1
                                total_unique_miseq[var][alt] = count

                for var in hiseq_var.keys():
                    for alt in hiseq_var[var]:
                        try:
                            if alt in miseq_var[var]:
                                pass
                            else:
                                sample_hiseq += 1
                                if var not in total_unique_hiseq:
                                    total_unique_hiseq[var] = {}
                                    total_unique_hiseq[var][alt] = 1
                                    count_unique_hiseq += 1
                                elif alt not in total_unique_hiseq[var]:
                                    total_unique_hiseq[var][alt] = 1
                                    count_unique_hiseq += 1
                                else:
                                    count = total_unique_hiseq[var][alt]
                                    count += 1
                                    total_unique_hiseq[var][alt] = count
                        except KeyError:
                            sample_hiseq += 1
                            if var not in total_unique_hiseq:
                                total_unique_hiseq[var] = {}
                                total_unique_hiseq[var][alt] = 1
                                count_unique_hiseq += 1
                            elif alt not in total_unique_hiseq[var]:
                                total_unique_hiseq[var][alt] = 1
                                count_unique_hiseq += 1
                            else:
                                count = total_unique_hiseq[var][alt]
                                count += 1
                                total_unique_hiseq[var][alt] = count
                sample_breakdown[v] = {'shared': sample_shared, 'miseq': sample_miseq, 'hiseq': sample_hiseq}
        else:
            miseq = [line.strip('\n') for line in open(fields[1], 'r').readlines()]
            miseq_var = generate_v_dict(miseq)
            hiseq = [line.strip('\n') for line in open(fields[2], 'r').readlines()]
            hiseq_var = generate_v_dict(hiseq)

            sample_shared = 0
            sample_miseq = 0
            sample_hiseq = 0

            for var in miseq_var.keys():
                for alt in miseq_var[var]:
                    try:
                        if alt in hiseq_var[var]:
                            sample_shared += 1
                            if var not in total_unique_shared:
                                total_unique_shared[var] = {}
                                total_unique_shared[var][alt] = 1
                                count_unique_shared += 1
                            elif alt not in total_unique_shared[var]:
                                total_unique_shared[var][alt] = 1
                                count_unique_shared += 1
                            else:
                                count = total_unique_shared[var][alt]
                                count += 1
                                total_unique_shared[var][alt] = count
                        else:
                            sample_miseq += 1
                            if var not in total_unique_miseq:
                                total_unique_miseq[var] = {}
                                total_unique_miseq[var][alt] = 1
                                count_unique_miseq += 1
                            elif alt not in total_unique_miseq[var]:
                                total_unique_miseq[var][alt] = 1
                                count_unique_miseq += 1
                            else:
                                count = total_unique_miseq[var][alt]
                                count += 1
                                total_unique_miseq[var][alt] = count

                    except KeyError:
                        sample_miseq += 1
                        if var not in total_unique_miseq:
                            total_unique_miseq[var] = {}
                            total_unique_miseq[var][alt] = 1
                            count_unique_miseq += 1
                        elif alt not in total_unique_miseq[var]:
                            total_unique_miseq[var][alt] = 1
                            count_unique_miseq += 1
                        else:
                            count = total_unique_miseq[var][alt]
                            count += 1
                            total_unique_miseq[var][alt] = count

            for var in hiseq_var.keys():
                for alt in hiseq_var[var]:
                    try:
                        if alt in miseq_var[var]:
                            pass
                        else:
                            sample_hiseq += 1
                            if var not in total_unique_hiseq:
                                total_unique_hiseq[var] = {}
                                total_unique_hiseq[var][alt] = 1
                                count_unique_hiseq += 1
                            elif alt not in total_unique_hiseq[var]:
                                total_unique_hiseq[var][alt] = 1
                                count_unique_hiseq += 1
                            else:
                                count = total_unique_hiseq[var][alt]
                                count += 1
                                total_unique_hiseq[var][alt] = count
                    except KeyError:
                        sample_hiseq += 1
                        if var not in total_unique_hiseq:
                            total_unique_hiseq[var] = {}
                            total_unique_hiseq[var][alt] = 1
                            count_unique_hiseq += 1
                        elif alt not in total_unique_hiseq[var]:
                            total_unique_hiseq[var][alt] = 1
                            count_unique_hiseq += 1
                        else:
                            count = total_unique_hiseq[var][alt]
                            count += 1
                            total_unique_hiseq[var][alt] = count

            sample_breakdown[s] = {'shared':sample_shared,'miseq':sample_miseq, 'hiseq':sample_hiseq}

    sample_breakdown['total'] = {'shared':count_unique_shared, 'hiseq':count_unique_hiseq, 'miseq':count_unique_miseq}

    f = open(args.out + '/unique_variants_summary.json', 'w')
    j = json.dumps(sample_breakdown, indent=4)
    print >> f, j
    f.close()

    f = open(args.out + '/unique_shared_variants.json', 'w')
    j = json.dumps(total_unique_shared, indent=4)
    print >> f, j
    f.close()

    f = open(args.out + '/unique_hiseq_variants.json', 'w')
    j = json.dumps(total_unique_hiseq, indent=4)
    print >> f, j
    f.close()

    f = open(args.out + '/unique_miseq_variants.json', 'w')
    j = json.dumps(total_unique_miseq, indent=4)
    print >> f, j
    f.close()

if __name__ == "__main__":
    variant_counts()













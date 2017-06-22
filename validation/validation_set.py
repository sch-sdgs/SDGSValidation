import argparse
import os
import subprocess
import glob
import json

def prepare_indel_vcf(vcf):
    """

    :param vcf:
    :return:
    """
    decomposed = vcf.replace('.vcf', '.decomposed.vcf')
    normalised = vcf.replace('.vcf', '.decomposed.normalised.vcf')
    indels_only = vcf.replace('.vcf', '.indels_only.vcf')

    try:
        command = '/results/Pipeline/program/vt/vt decompose ' + vcf + ' -o ' + decomposed
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        print(command)
        exit(1)

    try:
        command = '/results/Pipeline/program/vt/vt normalize ' + decomposed + ' -r ' + \
                  '/results/Pipeline/program/GATK_resource_bundle/ucsc.hg19.nohap.masked.fasta -o ' + normalised
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command:' + str(e.returncode))
        print(command)
        exit(1)

    try:
        command = '/results/Pipeline/program/bcftools-1.3.1/bcftools view -v "indels" ' + normalised + ' > ' + indels_only
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    print('vcf prepared for isec')
    return indels_only

def get_indel_summary(out_dir, vcf, folder, v1, cov, cov_file, final_stats):
    """

    :param out_dir:
    :param vcf:
    :param folder:
    :param v1:
    :param cov:
    :return:
    """
    fields = folder.split('_')
    pool = fields[0].replace('Pool', '')
    rd = fields[1]

    try:
        command = 'bgzip ' + vcf
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    zipped = vcf + '.gz'

    try:
        command = 'tabix ' + zipped
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        exit(1)

    stats = out_dir + folder + '/stats'
    command = '/results/Pipeline/program/bcftools-1.3.1/bcftools stats ' + v1 + ' ' + zipped + ' > ' + stats
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    command = 'grep "number of indels" ' + stats + ' | grep -v "\#" | cut -f1,2,4 | sed \'s/SN\t2/\tShared/g\' |  sed \'s/SN\t0/\tTruthOnly/g\'| sed \'s/SN\t1/\tPipelineOnly/g\' '
    try:
        out = subprocess.check_output(command, shell=True)
        lines = out.split('\n')
        for line in lines:
            if line == "":
                continue
            new_line = pool + '\t' + rd + line + '\n'
            final_stats.write(new_line)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    command = 'awk \' $10 != "" {print FILENAME"\t"$10}\' ' + out_dir + folder + '/0003.vcf | grep  ":" | cut -f1,3 -d":" | sed \'s/0\\/1://g\'|sed \'s/1\\/1://g\' >> ' + out_dir + 'depth_for_true_pos'
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    command = 'grep -v "#" ' + out_dir + folder + '/0001.vcf | cut -f6 >> ' + out_dir + 'quality_for_false_pos'
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    command = 'grep -v "#" ' + out_dir + folder + '/0003.vcf | cut -f6 >> ' + out_dir + 'quality_for_true_pos'
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print('Error executing command: ' + str(e.returncode))
        print(command)
        exit(1)

    f=open(cov, 'r')
    lines = [line.strip() for line in f.readlines()]
    f.close()

    for line in lines:
        if "median" in line:
            median = line.split(' ')[1].replace(',', '')
            print(median)
            cov_line = '\t'.join([pool,rd,median]) + '\n'
            cov_file.write(cov_line)


def main():
    """

    :return:
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', help='Location of output')
    parser.add_argument('-p', help='Sample name for reference genome as it appears in the VCF', default='S1613415-02')
    parser.add_argument('-b', help='Prefix associated with panel BED files (i.e. NGD_)')
    parser.add_argument('-bam', help='BAM file for run')
    parser.add_argument('-v', help='VCF to be compared')
    parser.add_argument('-rv', help='Path to reference genome vcf',
                        default='/results/Analysis/HiSeq_validation/giab/giab-NA12878/truth_small_variants.decomposed.normalised.vcf.gz')
    parser.add_argument('-rs', help='Identifier used in reference vcf', default='INTEGRATION')
    parser.add_argument('-indel', help='csv file listing all indel samples to be compared (v1,s1,v2,s2,bed,bam,cov,folder)')
    parser.add_argument('-mega', help='csv file listing all the mega samples to be compared (v1,s1,v2,s2,bed,bam,folder)')

    args = parser.parse_args()

    out_dir = args.o

    if not out_dir.endswith('/'):
        out_dir += '/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir + "GIAB"):
        os.makedirs(out_dir + "GIAB")

    #giab_comparison.giab_comp(out_dir + "GIAB/",args.p,args.b,args.bam,args.v,args.rv,args.rs)

    # f = open(args.indel, 'r')
    # files = [line.strip() for line in f.readlines()]
    # f.close()
    #
    # cov_file = open(out_dir + 'coverage', 'w')
    # cov_file.write('\t'.join(['pool', 'reads', 'median']) + '\n')
    # stats_file = open(out_dir + 'final_stats', 'w')
    #
    # for f in files:
    #     v1,s1,v2,s2,bed,bam,cov,folder = f.split(',')
    #     vcf = prepare_indel_vcf(v2)
    #     if not os.path.exists(out_dir + folder):
    #         os.makedirs(out_dir + folder)
    #
    #     command='python /home/bioinfo/Natalie/wc/sdgsValidation/validation/isec_vcf.py -o ' + out_dir + folder + \
    #             '/ -v1 ' + v1 + ' -v1_prep False -s1 ' + s1 + ' -v2 ' + vcf + ' -s2 ' + s2 + ' -b ' + bed + ' -bam ' + bam
    #
    #     try:
    #         subprocess.check_call(command, shell=True)
    #     except subprocess.CalledProcessError as e:
    #         print(command)
    #         print('Error executing command: ' + str(e.returncode))
    #         exit(1)
    #
    #     get_indel_summary(out_dir, vcf, folder, v1, cov, cov_file, stats_file)
    #
    # cov_file.close()
    # stats_file.close()
    #
    f = open(args.mega, 'r')
    files = [line.strip() for line in f.readlines()]
    f.close()

    for f in files:
        v1, s1, v2, s2, bed, bam, folder = f.split(',')
        if not os.path.exists(out_dir + folder):
            os.makedirs(out_dir + folder)

        command = 'python /home/bioinfo/Natalie/wc/sdgsValidation/validation/isec_vcf.py -o ' + out_dir + folder + \
                  '/ -v1 ' + v1 + ' -s1 ' + s1 + ' -v2 ' + v2 + ' -s2 ' + s2 + ' -b ' + bed + ' -bam ' + bam

        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            print(command)
            print('Error executing command: ' + str(e.returncode))
            exit(1)

    #indel summary
    # command = '/usr/bin/Rscript /home/bioinfo/Natalie/wc/sdgsValidation/validation/indel_summary.R ' + out_dir
    # try:
    #     subprocess.check_call(command, shell=True)
    # except subprocess.CalledProcessError as e:
    #     print(command)
    #     print('Error executing command: ' + str(e.returncode))
    #     exit(1)

    all_res = {}
    results = glob.glob(out_dir + '*/summary.json')
    for r in sorted(results):
        print(r)
        name = os.path.basename(os.path.dirname(r))
        print(name)
        f = open(r, 'r')
        d = json.load(f)
        f.close()
        out = d['result']
        all_res[name] = out

    f = open(out_dir + 'all_results.json', 'w')
    json.dump(all_res, f, indent=4)
    f.close()


if __name__ == "__main__":
    main()


from validation import SDGSvalidation
import argparse

def validate(baseline_file,baseline_bam,test_file,test_bam,sample,bed):
    log = SDGSvalidation.validate_variants(baseline_file,baseline_bam,test_file,test_bam,sample,bed)

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--baseline_file')
    parser.add_argument('--baseline_bam')
    parser.add_argument('--test_file')
    parser.add_argument('--test_bam')
    parser.add_argument('--sample')
    parser.add_argument('--bed',default=None)

    args = parser.parse_args()

    validate(args.baseline_file,args.baseline_bam,args.test_file,args.test_bam,args.sample,args.bed)

if __name__ == '__main__':
    main()
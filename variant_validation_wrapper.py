import os
import argparse
import glob
import json
from validation import misc_validation
import pandas

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f, 1):
            pass
    return i

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--worklist')

    old_version="v2"
    new_version="v3_1_1"

    args = parser.parse_args()

    samples = [x for x in os.listdir(args.worklist)]
    print samples
    result = {}

    for i in samples:
        if "S1607977" not in i:
            if i.startswith("S"):
                if os.path.isdir(args.worklist + "/" + i):
                    if new_version not in i:
                        sample=i.replace(old_version,"")
                        result[sample] = {}
    print result
    groups=[]
    for i in samples:
        if "S1607977" not in i:
            if i.startswith("S"):
                if os.path.isdir(args.worklist+"/"+i):
                    if old_version in i:
                        if new_version not in i:
                            sample=i.replace(old_version,"")
                            result[sample][old_version] = {}
                            bam=glob.glob(args.worklist+"/"+i+"/"+"*.bam")[0]
                            result[sample][old_version]["bam"] = bam
                            abbrev=''
                            log = glob.glob(args.worklist+"/"+i+"/"+"*_analysis_log.txt")
                            with open(log[0]) as f:
                                for l in f:
                                    if l.startswith("small-panel BED file"):
                                        line = l.rstrip().split(" ")
                                        bed = line[3]
                                        result[sample][old_version]["bed"]=bed
                                    if l.startswith("Disease group"):
                                        line = l.rstrip().split(" ")
                                        group = line [2]
                                        result[sample]["group"] = group
                                        if group not in groups:
                                            groups.append(group)
                                    if l.startswith("Abbreviation for small panel"):
                                        line=l.rstrip().split(" ")
                                        abbrev = line[4]
                            if len(groups) == 0:
                                result[sample]["group"] = 'PKD'
                                result[sample][old_version]["bed"] = '/results/Analysis/Ion_Torrent/BED_files/PKD.bed'
                            if not abbrev:
                                abbrev=''
                            variants = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/" + "*LessLQsPolys.txt")[0]
                            print variants

                            result[sample][old_version]["variants"] = variants

                            old_gaps = glob.glob(args.worklist+"/"+i+"/"+abbrev+"/*_coverage_curve_*.txt")

                            total_gaps_old = file_len(old_gaps[0]) - 1

                            min_file = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/*MinCoverage.txt")
                            min=int()
                            max=int()
                            if len(min_file) > 0:
                                with open(min_file[0]) as f:
                                    for l in f:
                                        min = l.rstrip()
                            else:
                                min ='0'
                            max_file = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/*MaxCoverage.txt")
                            if len(max_file)>0:
                                with open(max_file[0]) as f:
                                    for l in f:
                                        max = l.rstrip()
                            else:
                                max = '0'

                            # new_gaps = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/*_gaps_in_sequencing.txt")
                            # if len(new_gaps) > 0:
                            #     total_gaps_old = file_len(new_gaps[0]) - 1
                            # else:
                            #     total_gaps_old = "NA"
                            #
                            # summary = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/*_coverage_summary.txt")
                            # if len(summary) > 0:
                            #     with open(summary[0]) as json_file:
                            #         json_data = json.load(json_file)
                            #         max = json_data["max"]
                            #         min = json_data["min"]
                            # else:
                            #     max = "NA"
                            #     min = "NA"

                            result[sample][old_version]["gaps"] = total_gaps_old
                            result[sample][old_version]["min"] = int(min)
                            result[sample][old_version]["max"] = int(max)

                            # summary = glob.glob(args.worklist + "/" + i + "/" + abbrev + "/*_samtools_stats_all.json")
                            # if len(summary) > 0:
                            #     with open(summary[0]) as json_file:
                            #         json_data = json.load(json_file)
                            #         reads_mapped_and_paired = json_data["reads_mapped_and_paired"]

                            summary = glob.glob(args.worklist + "/" + i + "/*Alignment.txt")
                            if len(summary)>0:
                                with open(summary[0]) as f:
                                    for l in f:
                                        if "properly paired" in l:
                                            reads_mapped_and_paired=l.split(" ")[0]

                            else:
                                reads_mapped_and_paired=0
                            result[sample][old_version]["mapped_and_paired"] = reads_mapped_and_paired
                            summary = glob.glob(args.worklist + "/" + i + "/*_ReadsOffTarget.txt")
                            if len(summary)>0:
                                with open(summary[0]) as f:
                                    for l in f:
                                        result[sample][old_version]["off-target"] = l.rstrip()
                            else:
                                result[sample][old_version]["off-target"] = 0



                    else:
                        print i
                        sample, version = i.split("v")
                        result[sample][new_version] = {}
                        if "S" in i:
                            bam = glob.glob(args.worklist + "/" + i + "/" + "*.bam")[0]
                            result[sample][new_version]["bam"] = bam

                            log = glob.glob(args.worklist + "/" + i + "/" + "*_analysis_log.txt")
                            with open(log[0]) as f:
                                for l in f:
                                    if l.startswith("small-panel BED file"):
                                        line = l.rstrip().split(" ")
                                        bed = line[3]
                                        result[sample][new_version]["bed"] = bed
                            if "bed" not in result[sample][new_version]:
                                result[sample][new_version]["bed"] = '/results/Analysis/Ion_Torrent/BED_files/PKD.bed'
                            folder = [x for x in os.listdir(args.worklist + "/" + i)]
                            for analysis_name in folder:
                                if os.path.isdir(args.worklist + "/" + i + "/" + analysis_name):
                                    analysis_folder = os.path.basename(analysis_name)

                            variants = glob.glob(args.worklist + "/" + i + "/" + analysis_folder + "/" + "*LessLQsPolys.txt")[0]
                            result[sample][new_version]["variants"] = variants

                            new_gaps = glob.glob(args.worklist + "/" + i + "/" + analysis_folder + "/*_gaps_in_sequencing.txt")
                            if len(new_gaps) > 0:
                                total_gaps_new = file_len(new_gaps[0]) - 1
                            else:
                                total_gaps_new = "NA"

                            result[sample][new_version]["gaps"]=total_gaps_new

                            summary = glob.glob(args.worklist + "/" + i + "/" + analysis_folder + "/*_coverage_summary.txt")
                            print summary
                            if len(summary) > 0:
                                with open(summary[0]) as json_file:
                                    json_data = json.load(json_file)
                                    max = json_data["max"]
                                    min = json_data["min"]
                            else:
                                max = "NA"
                                min = "NA"

                            result[sample][new_version]["min"] = min
                            result[sample][new_version]["max"] = max

                            summary = glob.glob(args.worklist + "/" + i + "/*Alignment.txt")
                            if len(summary)>0:
                                with open(summary[0]) as f:
                                    for l in f:
                                        if "properly paired" in l:
                                            reads_mapped_and_paired = l.split(" ")[0]
                            else:
                                reads_mapped_and_paired = 0

                            result[sample][new_version]["mapped_and_paired"] = reads_mapped_and_paired

                            summary = glob.glob(args.worklist + "/" + i + "/*_ReadsOffTarget.txt")
                            if len(summary) > 0:
                                with open(summary[0]) as f:
                                    for l in f:
                                        result[sample][new_version]["off-target"] = l.rstrip()
                            else:
                                result[sample][new_version]["off-target"] = 0

    print json.dumps(result,indent=4)
    for sample in result:
        print sample
        old=result[sample][old_version]
        new=result[sample][new_version]
        print old["bam"]
        print new["bam"]
        variant_validation = misc_validation.validate_variants(old["variants"], old["bam"], new["variants"], new["bam"], sample,"/results/Pipeline/masterBED_Backup/NGD_ataxia_150bp.bed")


        variants = variant_validation["missing_variants"].to_dict(orient="records")
        missing_variants = []
        for variant in variants:
            chrom = variant["chrom"]
            pos = variant["pos"]
            ref = variant["ref"]
            alt = variant["alt"]
            qual = variant["qual"]
            data = ''
            for j in variant["info"].split(";"):
                key,value=j.split("=")
                if key == "dataset":
                    data = value
            var = [chrom,pos,ref,alt,qual,data]
            missing_variants.append(":".join(var))

        variants = variant_validation["additional_variants"].to_dict(orient="records")
        additional_variants = []
        for variant in variants:
            chrom = variant["chrom"]
            pos = variant["pos"]
            ref = variant["ref"]
            alt = variant["alt"]
            qual = variant["qual"]
            data = ''
            for j in variant["info"].split(";"):
                key, value = j.split("=")
                if key == "dataset":
                    data = value
            var = [chrom, pos, ref, alt, qual, data]
            additional_variants.append(":".join(var))
        new["variants_missing_actual"] = missing_variants
        old["variants_missing_actual"] = additional_variants
        data = variant_validation["stats"].to_dict()
        shared = data["Variants Shared"][0]
        missing_old,x = data["Pass Variants Missing (All)"][0].split(" ")
        missing_new,x = data["Pass Variants Missing (All)"][1].split(" ")
        new["variants_missing"] = missing_new
        old["variants_missing"] = missing_old
        old["variants_shared"] = shared
        new["variants_shared"] = shared
        if new["max"] != "NA":
            max_diff=new["max"] - old["max"]
        else:
            max_diff="NA"
        if new["min"] != "NA":
            min_diff = new["min"] - old["min"]
        else:
            min_diff = "NA"
        if new["gaps"] != "NA":
            gaps_diff=new["gaps"] - old["gaps"]
        else:
            gaps_diff = "NA"

        result[sample]["max_diff"] = max_diff
        result[sample]["min_diff"] = min_diff
        result[sample]["gaps_diff"] = gaps_diff

        #do bedtools intersect to create a file that gives the new sites
        if gaps_diff > 0:
            pass
    if len(groups) == 0:
        groups.append('PKD')
    for g in groups:


        header=["sample","disease_group",old_version+"_bed",new_version+"_bed","bed_status",old_version+"_mapped",new_version+"_mapped","expected_mapped",old_version+"_min",new_version+"_min","min_diff","min_status",old_version+"_max",new_version+"_max","max_diff",old_version+"_gaps",new_version+"_gaps","gaps_diff","gap_status","shared_variants","missing_in_"+old_version,"missing_in_"+new_version,"variant_status","var_missing_in_"+old_version,"var_missing_in_"+new_version+"\n"]

        file = open(args.worklist+"/"+g+".csv",'w')

        file.write(",".join(header))


        for sample in result:
            print json.dumps(result,indent=4)
            print sample
            group=result[sample]["group"]
            if group == g:
                v3_bed=result[sample][old_version]["bed"]
                v3_1_bed = result[sample][new_version]["bed"]
                v3_missing=result[sample][old_version]["variants_missing"]
                v3_1_missing=result[sample][new_version]["variants_missing"]
                shared=result[sample][old_version]["variants_shared"]
                variant_status = "CHECK"
                if v3_1_missing == "0":
                    if v3_missing == "0":
                        variant_status = "OK"
                bed_status = "OK"
                if v3_bed != v3_1_bed:
                    bed_status="BED_CHANGED"
                v3_min=result[sample][old_version]["min"]
                v3_1_min=result[sample][new_version]["min"]
                min_status = "OK"
                if v3_min >= 30:
                    if v3_1_min < 30:
                        min_status = "CHECK"
                v3_max = result[sample][old_version]["max"]
                v3_1_max = result[sample][new_version]["max"]
                v3_gaps = result[sample][old_version]["gaps"]
                v3_1_gaps = result[sample][new_version]["gaps"]
                v3_mapped = result[sample][old_version]["mapped_and_paired"]
                v3_1_mapped = result[sample][new_version]["mapped_and_paired"]
                max_diff = result[sample]["max_diff"]
                min_diff = result[sample]["min_diff"]
                gaps_diff = result[sample]["gaps_diff"]
                gap_status = "OK"
                if gaps_diff > 0:
                    gap_status = "CHECK"
                var_missing_new = "|".join(result[sample][new_version]["variants_missing_actual"])
                var_missing_old = "|".join(result[sample][old_version]["variants_missing_actual"])

                approx_expected_reads=int(v3_mapped)-int(result[sample][old_version]["off-target"])

                out = [sample,group,v3_bed,v3_1_bed,bed_status,v3_mapped,v3_1_mapped,approx_expected_reads,v3_min,v3_1_min,min_diff,min_status,v3_max,v3_1_max,max_diff,v3_gaps,v3_1_gaps,gaps_diff,gap_status,shared,v3_missing,v3_1_missing,variant_status,var_missing_old,var_missing_new+"\n"]
                file.write(",".join(str(x) for x in out))
        file.close()
    #print json.dumps(result, indent=4)

if __name__ == '__main__':
    main()
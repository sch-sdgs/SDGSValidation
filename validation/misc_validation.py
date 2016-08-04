from producers import VariantProducers
import json
import logging
import pandas as pd
import sqlite3
import time
# from weasyprint import HTML
from pybedtools import BedTool

logger = logging.getLogger(__name__)

logging.basicConfig(level=logging.ERROR)


def determine_file_type(file):
    f = open(file)
    head = f.readline()
    if head.startswith('Variant'):
        producer = VariantProducers.ProcessedVariantsV1
        return producer
    if "cons_preferred_transcripts" in head:
        producer = VariantProducers.ProcessedVariantsV3
        return producer
    if head.startswith('chrom'):
        producer = VariantProducers.ProcessedVariantsV2
        return producer
    if head.startswith('##fileformat=VCFv4.1'):
        producer = VariantProducers.GVCF
        return producer
    if head.startswith('##fileformat=VCFv4.2'):
        producer = VariantProducers.GVCF
        return producer



def check_match(baseline, test):
    if type(baseline) == str or type(baseline) == int:
        if baseline == test:
            return True
        else:
            return False
    elif type(baseline) == list:
        for i in test:
            if i in baseline:
                return True
            else:
                return False


def append_variants_to_table(variant, table, additional_field=None, additional_value=None):
    for_table = str(variant).rstrip("\n").split("\t")
    if additional_field is not None:
        for i, field in enumerate(additional_field):
            if field not in table:
                table[field] = []
            table[field].append(additional_value[i])

    table["chrom"].append(for_table[0])
    table["pos"].append(for_table[1])
    table["id"].append(for_table[2])
    table["ref"].append(for_table[3])
    table["alt"].append(for_table[4])
    table["qual"].append(for_table[5])
    table["filter"].append(for_table[6])
    table["info"].append(for_table[7])
    table["key"].append(for_table[8])
    table["genotype"].append(for_table[9])

    return table


def test_mismatch_list(test_list, baseline_list, parameter, baseline_variant, test_variant):
    """
    tests if two lists match each other

    :param test_list:
    :param baseline_list:
    :param parameter:
    :param baseline_variant:
    :param test_variant:
    :return:
    """
    if len(test_list) + len(baseline_list) != 0:
        if not any(map(lambda v: v in sorted(test_list), sorted(baseline_list))):
            logger.info(parameter + ": " + str(baseline_list) + " vs " + str(test_list))
            logger.info("   BASELINE_VARIANT: " + str(baseline_variant))
            logger.info("   TEST_VARIANT: " + str(test_variant))
            return False


def make_igv_link(baseline_bam, test_bam, chrom, pos):
    link = '<a href = "http://localhost:60151/goto?sort=base&locus=' + chrom + ':' + str(
        pos + 1) + '"><button type="button" class="btn btn-default"><span class="glyphicon glyphicon-search"></span></button></a>'
    return link


def test_mismatch_list(test_list, baseline_list, parameter, baseline_variant, test_variant):
    """
    tests if two lists match each other

    :param test_list:
    :param baseline_list:
    :param parameter:
    :param baseline_variant:
    :param test_variant:
    :return:
    """
    if len(test_list) + len(baseline_list) != 0:
        if not any(map(lambda v: v in sorted(test_list), sorted(baseline_list))):
            logger.info(parameter + ": " + str(baseline_list) + " vs " + str(test_list))
            logger.info("   BASELINE_VARIANT: " + str(baseline_variant))
            logger.info("   TEST_VARIANT: " + str(test_variant))
            return False

def select_group(id,type):

    select = "<p><input size=\"50\" class=\"form-control inputwide\" list=\"reasonList\" type=\"text\" name=\"variant_"+id+"-"+type+"\"> \
    <datalist id=\"reasonList\"> \
    <option value=\"Strand Bias\"> \
    <option value=\"Coverage\"> \
    <option value=\"Artefact\"> \
    <option value=\"Low Quality\"> \
    </datalist></p>"

    return select

def validate_variants(baseline_file, baseline_bam, test_file, test_bam, sample_id, bed_file):
    producer = determine_file_type(baseline_file)
    variant_producer_baseline = producer(baseline_file, sample_id)

    producer = determine_file_type(test_file)
    variant_producer_test = producer(test_file, sample_id)

    logger.info("Baseline File: " + baseline_file)
    logger.info("Test File: " + test_file)

    baseline_variants_by_id = {}
    test_variants_by_id = {}
    baseline_additional = {}
    test_additional = {}

    assert cmp(baseline_variants_by_id, test_variants_by_id) == 0

    for i in variant_producer_baseline.get_variants(bed=bed_file):
        variant = i.toJsonDict()
        baseline_variants_by_id[variant["id"]]=i


    for i in variant_producer_test.get_variants(bed=bed_file):
        variant = i.toJsonDict()
        test_variants_by_id[variant["id"]] = i


    logger.info("Baseline Variants Count: " + str(len(baseline_variants_by_id)))
    logger.info("Test Variants Count: " + str(len(test_variants_by_id)))

    shared_variants = [k for k in baseline_variants_by_id if k in test_variants_by_id]
    missing = [k for k in baseline_variants_by_id if k not in test_variants_by_id]

    for m in missing:
        logger.info("MISSING_VARIANT: " + str(baseline_variants_by_id[m]))

    # TODO: for vcf filter on lowquality flag in filter column

    fail_reasons = {"MISMATCH": 0}

    mismatch_variants = {"chrom": [], "pos": [], "id": [], "ref": [], "alt": [], "qual": [], "filter": [], "info": [],
                         "key": [], "genotype": []}
    additional_variants = {"chrom": [], "pos": [], "id": [], "ref": [], "alt": [], "qual": [], "filter": [], "info": [],
                           "key": [], "genotype": []}
    missing_variants = {"chrom": [], "pos": [], "id": [], "ref": [], "alt": [], "qual": [], "filter": [], "info": [],
                        "key": [], "genotype": []}


    for i in missing:
        variant = baseline_variants_by_id[i].toJsonDict()
        link = make_igv_link(baseline_bam, test_bam, variant["referenceName"], variant["start"])
        missing_variants = append_variants_to_table(baseline_variants_by_id[i], missing_variants,
                                                    additional_field=["reason", "bam"], additional_value=[select_group(i,"missing"), link])

    for i in shared_variants:
        baseline_variant = baseline_variants_by_id[i].toJsonDict()
        test_variant = test_variants_by_id[i].toJsonDict()

        baseline_alt = baseline_variant["alternateBases"]
        test_alt = test_variant["alternateBases"]

        baseline_genotype = baseline_variant["calls"][0]["genotype"]
        test_genotype = test_variant["calls"][0]["genotype"]

        baseline_info = baseline_variant["info"]
        test_info = test_variant["info"]

        if "transcripts" in baseline_info:
            baseline_tx = baseline_info["transcripts"]
        else:
            baseline_tx = []

        if "transcripts" in test_info:
            test_tx = test_info["transcripts"]
        else:
            test_tx = []

        baseline_rsid = baseline_variant["names"]
        test_rsid = test_variant["names"]

        genotype_check = test_mismatch_list(test_genotype, baseline_genotype, "GENOTYPE_MISMATCH",
                                            baseline_variants_by_id[i], test_variants_by_id[i])
        alt_check = test_mismatch_list(test_alt, baseline_alt, "ALT_MISMATCH", baseline_variants_by_id[i],
                                       test_variants_by_id[i])
        tx_check = test_mismatch_list(test_tx, baseline_tx, "TX_MISMATCH", baseline_variants_by_id[i],
                                      test_variants_by_id[i])
        rsid_check = test_mismatch_list(test_rsid, baseline_rsid, "RSID_MISMATCH", baseline_variants_by_id[i],
                                        test_variants_by_id[i])

        link = make_igv_link(baseline_bam, test_bam, baseline_variant["referenceName"], baseline_variant["start"])

        if genotype_check == False or alt_check == False or tx_check == False or rsid_check == False:
            mismatch_variants = append_variants_to_table(baseline_variants_by_id[i], mismatch_variants,
                                                         additional_field=["scope", "mismatch", "reason", "bam"],
                                                         additional_value=["BASELINE", "", "", ""])
            mismatch_variants = append_variants_to_table(test_variants_by_id[i], mismatch_variants,
                                                         additional_field=["scope", "mismatch", "reason", "bam"],
                                                         additional_value=["TESTING", "",select_group(i,"mismatch"),link])
            fail_reasons["MISMATCH"] += 1

    additional_calls = [k for k in test_variants_by_id if k not in baseline_variants_by_id]
    for i in additional_calls:
        variant = test_variants_by_id[i]
        link = make_igv_link(baseline_bam, test_bam, variant.toJsonDict()["referenceName"],
                             variant.toJsonDict()["start"])
        additional_variants = append_variants_to_table(variant, additional_variants, additional_field=["reason", "bam"],
                                                       additional_value=[select_group(i,"additional"), link])
        logger.info("ADDITIONAL_VARIANT: " + str(variant))


    pd.set_option('display.max_colwidth', -1)


    additional_variants_pd = pd.DataFrame(additional_variants)
    missing_variants_pd = pd.DataFrame(missing_variants)


    mismatch_variants_pd = pd.DataFrame(mismatch_variants)

    if len(additional_variants['chrom']) > 0:
        additional_variants_pd = additional_variants_pd[
            ['bam', 'reason', 'chrom', 'pos', 'ref', 'alt', 'filter', 'key', 'genotype', 'id', 'qual', 'info']]

    if len(mismatch_variants['chrom']) > 0:
        mismatch_variants_pd = mismatch_variants_pd[
            ['bam', 'reason', 'scope', 'chrom', 'pos', 'ref', 'alt', 'filter', 'key', 'genotype', 'id', 'qual', 'info']]
        mismatch_status = 'FAIL'
    else:
        mismatch_status = 'PASS'
    if len(missing_variants['chrom']) > 0:
        missing_variants_pd = missing_variants_pd[
            ['bam', 'reason', 'chrom', 'pos', 'ref', 'alt', 'filter', 'key', 'genotype', 'id', 'qual', 'info']]
    # print missing_variants_pd
    # missing_variants_pd = missing_variants_pd.loc[missing_variants_pd['filter'] != 'LowQual']
    missing_variants_pass_count = len(missing_variants_pd.index)
    # mismatch_variants_pd = mismatch_variants_pd.loc[mismatch_variants_pd['filter'] == "PASS"]
    # additional_variants_pd = additional_variants_pd.loc[additional_variants_pd['filter'] != 'LowQual']
    additional_variants_pass_count = len(additional_variants_pd.index)

    stats = pd.DataFrame({
        'Type': ["baseline", "test"],
        'File': [baseline_file, test_file],
        'Total Variants': [len(baseline_variants_by_id), len(test_variants_by_id)],
        'Pass Variants Missing (All)': [str(additional_variants_pass_count) + " (" + str(len(additional_calls)) + ")",
                                        str(missing_variants_pass_count) + " (" + str(len(missing)) + ")"],
        'Variants Shared': [len(shared_variants), len(shared_variants)],
    'Variants Fail Mismatch': ['-', fail_reasons["MISMATCH"]]})


    template_vars = {"title": "Validation Report: " + sample_id,
                 "stats": stats,
                 "bam_link": "http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
                 "missing_variants": missing_variants_pd,
                 "mismatch_variants": mismatch_variants_pd,
                 "mismatch_status": mismatch_status,
                 "additional_variants": additional_variants_pd
                 }

    return template_vars


def validate_stats():
    pass
    # total reads
    # percentage dups
    # library size
    # median read length
    # percent mapped
    # correctly orientated
    # percentage reads in broad panel region
    # median insert size
    # percentage of small panel region with depth at least 30x

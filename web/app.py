from flask import Flask, render_template, Markup, request, url_for
from validation import SDGSvalidation
import os

app = Flask(__name__)


@app.route('/')
def form():
    in_progress = os.listdir("in_progress")
    complete = os.listdir("complete_validations")
    return render_template('form_submit.html',complete=complete,in_progress=in_progress)

@app.route('/hello/', methods=['POST'])
def hello():
    baseline_file=request.form['baseline_file']
    baseline_bam=request.form['baseline_bam']
    test_file = request.form['test_file']
    test_bam = request.form['test_bam']
    sample = request.form['sample']
    bed = request.form['bed']
    result = SDGSvalidation.validate_variants(baseline_file, baseline_bam, test_file, test_bam, sample, bed)
    return render_template('validation_report_variants.html',
                           title="Validation report",
                           stats=Markup(result["stats"]),
                           bam_link="http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
                           missing_variants=Markup(result["missing_variants"]),
                           mismatch_variants=Markup(result["mismatch_variants"]),
                           additional_variants=Markup(result["additional_variants"]))


# @app.route("/")
# def main():
#     baseline_file='/results/Analysis/HiSeq/2015/1504850/S0902255-07/Ataxia/1504850-S0902255-07_ConsensusVariants_Ataxia_LessLQs.txt'
#     baseline_bam='/results/Analysis/HiSeq/2015/1504850/S0902255-07/1504850-S0902255-07_Aligned_Sorted_PCRDuped_IndelsRealigned.bam'
#     test_file='/home/bioinfo/mparker/work/20160613/pipeline/1504850-S0902255-07_Ataxia_variants_with_downsampling.txt'
#     test_bam='/home/bioinfo/mparker/work/20160613/pipeline/1504850-S0902255-07_Aligned_Sorted_PCRDuped_IndelsRealigned.bam'
#     sample='TEST'
#     bed='/results/Pipeline/masterBED_Backup/NGD_ataxia_150bp.bed'
#     result = SDGSvalidation.validate_variants(baseline_file, baseline_bam, test_file, test_bam, sample, bed)
#     return render_template('validation_report_variants.html',
#                            title="Validation report",
#                            stats=Markup(result["stats"]),
#                            bam_link="http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
#                            missing_variants=Markup(result["missing_variants"]),
#                            mismatch_variants=Markup(result["mismatch_variants"]),
#                            additional_variants=Markup(result["additional_variants"]))

#
# template_vars = {"title": "Validation report: " + sample_id,
#                  "stats": stats.to_html(index=False, justify="left", classes=["table table-striped"]),
#                  "bam_link": "http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
#                  "missing_variants": missing_variants_pd.to_html(index=False, justify="left", escape=False,
#                                                                  classes=["table table-striped"]),
#                  "mismatch_variants": mismatch_variants_pd.to_html(index=False, justify="left", escape=False,
#                                                                    classes=["table table-striped"]),
#                  "mismatch_status": mismatch_status,
#                  "additional_variants": additional_variants_pd.to_html(index=False, justify="left", escape=False,
#                                                                        classes=["table table-striped"])}

if __name__ == "__main__":
    # app.run(host= '10.182.131.21')
    app.run()
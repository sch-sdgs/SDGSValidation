from flask import Flask, render_template, Markup, request, url_for
from validation import SDGSvalidation
import os

app = Flask(__name__)


@app.route('/')
def form():
    in_progress = os.listdir("in_progress")
    complete = os.listdir("complete_validations")
    return render_template('validation_home.html', complete=complete, in_progress=in_progress)


@app.route('/variant_validation/', methods=['POST'])
def variant_validation():
    baseline_file = request.form['baseline_file']
    baseline_bam = request.form['baseline_bam']
    test_file = request.form['test_file']
    test_bam = request.form['test_bam']
    sample = request.form['sample']
    bed = request.form['bed']
    result = SDGSvalidation.validate_variants(baseline_file, baseline_bam, test_file, test_bam, sample, bed)
    return render_template('validation_report_variants.html',
                           title="Validation Report",
                           stats=Markup(result["stats"]),
                           bam_link="http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
                           missing_variants=Markup(result["missing_variants"]),
                           mismatch_variants=Markup(result["mismatch_variants"]),
                           additional_variants=Markup(result["additional_variants"]))

if __name__ == "__main__":
    app.run(host= '10.182.131.21', port='5001')
    #app.run()

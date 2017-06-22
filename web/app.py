from flask import Flask, render_template, Markup, request, url_for
import sys
from validation import misc_validation
import os
import json
from flask.ext.wtf import Form
from wtforms.fields import TextField, SubmitField
import pprint
app = Flask(__name__)


@app.route('/')
def form():
    # paths=dict()
    # d = {}
    # miseq_runs = os.listdir('/results/Analysis/MiSeq/Runs')
    # for i in miseq_runs:
    #     if i.startswith('2016'):
    #         paths[i]={}
    #         workflows = os.listdir('/results/Analysis/MiSeq/Runs/'+i)
    #         for w in workflows:
    #             paths[i][w]={}
    #             samples = os.listdir('/results/Analysis/MiSeq/Runs/'+i+'/'+w)
    #             for s in samples:
    #                 if s.startswith('S'):
    #                     paths[i][w][s]=[]
    #                     path = '/results/Analysis/MiSeq/Runs/' + i + '/' + w + '/' + s
    #                     if os.path.isdir(path):
    #                         files = os.listdir(path)
    #                         for f in files:
    #                             path = '/results/Analysis/MiSeq/Runs/' + i + '/' + w + '/' + s + '/' + f
    #                             if os.path.isdir(path):
    #                                 d[f]=[]
    #                                 files_last = os.listdir(path)
    #                                 for l in files_last:
    #                                     d[f].append(l)
    #                                 paths[i][w][s].append(d)
    #                             else:
    #                                 paths[i][w][s].append(f)


    in_progress = os.listdir("in_progress")
    complete = os.listdir("complete_validations")
    return render_template('validation_home.html', complete=complete, in_progress=in_progress)

@app.route('/save_form/', methods=['POST'])
def save_form():
    if 'complete' in request.form:
        iscomplete = request.form['complete']
    else:
        iscomplete = False


    if 'status' in request.form:
        ispass = request.form['status']
    else:
        ispass = False

    for i in request.form:
        if i.startswith('variant'):
            print i
            print request.form[i]

    if iscomplete:
        return render_template('validation_complete.html', ispass=ispass)

@app.route('/variant_validation/', methods=['POST'])
def variant_validation():
    for i in request.form:
        print i
    baseline_file = request.form['baseline_file']
    baseline_bam = request.form['baseline_bam']
    test_file = request.form['test_file']
    test_bam = request.form['test_bam']
    sample = request.form['sample']
    bed = request.form['bed']
    result = misc_validation.validate_variants(baseline_file, baseline_bam, test_file, test_bam, sample, bed)

    return render_template('validation_report_variants.html',
                           title="Validation Report",
                           sample=sample,
                           stats=Markup(result["stats"].to_html(index=False, justify="left",
                                            classes=["table table-bordered table-striped"])),
                           bam_link="http://localhost:60151/load?file=http://10.182.131.21" + baseline_bam + ",http://10.182.131.21" + test_bam + "&merge=false",
                           missing_variants=Markup(result["missing_variants"].to_html(index=False, justify="left", escape=False,
                                                                     classes=["table table-bordered table-striped"])),
                           mismatch_variants=Markup(result["mismatch_variants"].to_html(index=False, justify="left", escape=False,
                                                                       classes=["table table-bordered table-striped"])),
                           additional_variants=Markup(result["additional_variants"].to_html(index=False, justify="left", escape=False,
                                                                           classes=[
                                                                               "table table-bordered table-striped"])))



@app.route('/giab/view', methods=["GET", "POST"])
def giab_view():
    if request.method == "GET":
        return render_template('view_giab.html')

@app.route('/giab/results', methods=["GET", "POST"])
def giab_results():
    print('requested')
    print(request.method)
    if request.method == "GET":
        print('requested')
        print(request.method)
        filepath = request.args.get("filepath")
        panel = request.args.get("panel")
        print(panel)
        print(filepath)
        f = open(filepath, 'r')
        j = json.load(f)
        f.close()
        return render_template('giab_results.html', panel=j[panel], panelname=panel, filepath=filepath)

@app.route('/giab/summary', methods=["GET", "POST"])
def giab_summary():
    print('requested')
    print(request.method)
    if request.method == "POST":
        print('requested')
        print(request.method)
        filepath = request.form["filepath"]
        print(filepath)
        f = open(filepath, 'r')
        j = json.load(f)
        f.close()
        return render_template('giab_summary.html', panels=j, order=sorted(j.keys()), filepath=filepath)

if __name__ == "__main__":
    app.run(host= '10.182.131.21',threaded=True,port=5010)
    app.secret_key = 'development key'





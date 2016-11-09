from __future__ import print_function
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
import json
import time
import os
import numpy as np


env = Environment(loader=FileSystemLoader('.'))



template = env.get_template("templates/standard_validation.html")

template_vars = {"title": "General: Protocol Validation Form SOP (New Test)",
    "section" : "Bioinformatics"}

html_out = template.render(template_vars)
HTML(string=html_out).write_pdf("/home/bioinfo/test.pdf",stylesheets=["static/simple_report.css"])
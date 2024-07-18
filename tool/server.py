import os
import subprocess
import random
import json
from Bio.Seq import Seq
from flask_wtf import FlaskForm, CSRFProtect, file
from werkzeug.utils import secure_filename
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired
from flask import Flask, render_template, current_app
from flask_wtf.file import FileField
from wtforms import SubmitField
from flask import Flask, render_template, request, redirect, url_for, abort
import re
from server_utils import *


# Initialize the Flask app
app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['SECRET_KEY'] = 'taltal'
csrf = CSRFProtect(app)  # Initialize CSRF protection
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024
app.config['UPLOAD_EXTENSIONS'] = ['.fasta']


class InputForm(FlaskForm):
    email = StringField("Email", validators=[DataRequired()], render_kw={"id": "email"})
    gene = StringField("Input gene", validators=[DataRequired()], render_kw={"id": "gene"})
    user_trigger = BooleanField("Got a known trigger?", render_kw={"id": "user_trigger"})
    trigger = StringField("Input trigger", render_kw={"id": "trigger"})
    rna = StringField("Input mRNA", render_kw={"id": "rna"})
    cell_type = SelectField("Organism Type", choices=[('Prokaryote', 'Prokaryote'), ('Eukaryote', 'Eukaryote')], render_kw={"id": "cell_type"})
    file = FileField('File', render_kw={"id": "file"})
    submit = SubmitField("Submit")



# Home page route
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')



@app.route('/form', methods=['GET', 'POST'])
def user_data_getter():
    email = None
    gene = None
    trigger = None
    rna = None
    cell_type = None
    input_form = InputForm()

    if input_form.validate_on_submit():
        trigger = input_form.trigger.data.upper()
        gene = input_form.gene.data.upper()
        rna = input_form.rna.data.upper()

        # Get the data from the form
        email = input_form.email.data
        gene = input_form.gene.data
        trigger = input_form.trigger.data
        mrna = input_form.rna.data
        cell_type = input_form.cell_type.data
        uploaded_file = input_form.file.data

        # TODO: change to process the file in run_scoring.py not in the user getter function?
        try:
            data_dict = generic_process_file(uploaded_file)
        except Exception as e:
            raise f"Error processing file: {e}"

        # assign and send the data to run_scoring.py in a subprocess
        s_email = str(email) if email else "EMPTY"
        s_gene = str(gene) if gene else "EMPTY"
        s_trigger = str(trigger) if trigger else "EMPTY"
        s_mrna = str(mrna) if mrna else "EMPTY"
        s_cell_type = str(cell_type) if cell_type else "EMPTY"
        s_file_dict = json.dumps(data_dict) if data_dict else "EMPTY"

        subprocess.Popen(["python3", "-m", "debugpy", "--listen", str(random.randint(10000, 50000)), "run_scoring.py",
                          s_gene, s_trigger, s_mrna, s_cell_type, s_email, s_file_dict])

        # Clear the form
        input_form.gene.data = ''
        input_form.trigger.data = ''
        input_form.rna.data = ''
        input_form.email.data = ''
        input_form.cell_type.data = ''
        input_form.file.data = ''

    return render_template('form.html', input_form=input_form, gene=gene, trigger=trigger, rna=rna, cell_type=cell_type, email=email, file=file)


# Error Handling##############################################################
# Invalid URL
@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404

# Internal Server Error
@app.errorhandler(500)
def internal_error(e):
    return render_template("500.html"), 500


if __name__ == '__main__':

    app.run(debug=True, port=3000)

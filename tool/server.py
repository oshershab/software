import os
import subprocess
import random
import json
from Bio.Seq import Seq
from wtforms import validators
from flask_wtf import FlaskForm, CSRFProtect
from werkzeug.utils import secure_filename
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired, Email, Length
from flask import Flask, render_template, request, redirect, url_for, flash
from flask_wtf.file import FileField
from server_utils import generic_process_file
from wtforms import StringField, TextAreaField, SubmitField


# Initialize the Flask app
app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['SECRET_KEY'] = 'taltal'
csrf = CSRFProtect(app)  # Initialize CSRF protection
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024
app.config['UPLOAD_EXTENSIONS'] = ['.fasta']

class InputForm(FlaskForm):
    # todo: add validators to email?
    email = StringField("Email", [DataRequired()], render_kw={"id": "email"})
    gene = StringField("Input Gene",  render_kw={"id": "gene"})
    user_trigger = BooleanField("Got a known trigger?",  render_kw={"id": "user_trigger"})
    trigger = StringField("Input Trigger", render_kw={"id": "trigger"})
    reporter_gene = StringField("Reporter Gene", [DataRequired()], render_kw={"id": "reporter_gene"})
    cell_type = SelectField("Organism Type", choices=[('Prokaryote', 'Prokaryote'), ('Eukaryote', 'Eukaryote')], render_kw={"id": "cell_type"})
    file = FileField('File', render_kw={"id": "file"})
    submit = SubmitField("Submit")

# Home page route
@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/form', methods=['GET', 'POST'])
def user_data_getter():
    print("In form")

    # TODO: validation before submission
    user_trigger_bool = False
    reporter_gene = None
    cell_type = None
    trigger = None
    email = None
    file = None
    gene = None




    input_form = InputForm()
    bool = True
    if bool:
        print("Form submitted")
        try:
            # Get the data from the form
            email = input_form.email.data
            if input_form.user_trigger.data:
                trigger = input_form.trigger.data.upper()
            else:
                gene = input_form.gene.data.upper()
            reporter_gene = input_form.reporter_gene.data.upper()
            cell_type = input_form.cell_type.data
            user_trigger_bool = input_form.user_trigger.data
            uploaded_file = input_form.file.data

            # Process the file
            data_dict = generic_process_file(uploaded_file)
            # Assign and send the data to preprocess.py in a subprocess
            s_email = str(email) if email else "EMPTY"
            s_gene = str(gene) if gene else "EMPTY"
            s_trigger = str(trigger) if trigger else "EMPTY"
            s_reporter = str(reporter_gene) if reporter_gene else "EMPTY"
            s_cell_type = str(cell_type) if cell_type else "EMPTY"
            s_file_dict = json.dumps(data_dict) if data_dict else "EMPTY"
            subprocess.Popen(["python3", "-m", "debugpy", "--listen", str(random.randint(10000, 50000)), "preprocess.py",
                              s_gene, s_trigger, s_reporter, s_cell_type, s_email, s_file_dict])

            # Clear the form
            input_form.gene.data = ''
            input_form.trigger.data = ''
            input_form.reporter_gene.data = ''
            input_form.email.data = ''
            input_form.cell_type.data = ''
            input_form.file.data = ''
            flash('Form submitted successfully. Job accepted.')

        except Exception as e:
            flash(f"Error processing form: {e}", "danger")

    return render_template('form.html', input_form=input_form)

# Error Handling
@app.errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404

@app.errorhandler(500)
def internal_error(e):
    return render_template("500.html"), 500

if __name__ == '__main__' :
    app.run(debug=True, port=5000)

import subprocess
from random import random
from flask_wtf import FlaskForm, CSRFProtect
from wtforms import StringField, BooleanField, SubmitField, SelectField
from wtforms.validators import DataRequired
from flask import Flask, render_template


app = Flask(__name__)
app.config['SECRET_KEY'] = 'taltal'
csrf = CSRFProtect(app)  # Initialize CSRF protection



# Home page route
@app.route('/')
def index():
    return render_template('index.html')



# Create a form class for submitting data
class InputForm(FlaskForm):
    email = StringField("Email", validators=[DataRequired()], render_kw={"id": "email"})
    gene = StringField("Input gene", validators=[DataRequired()], render_kw={"id": "gene"})
    user_trigger = BooleanField("Got a known trigger?", render_kw={"id": "user_trigger"})
    trigger = StringField("Input trigger", render_kw={"id": "trigger"})
    rna = StringField("Input mRNA", render_kw={"id": "rna"})
    cell_type = SelectField("Organism Type", choices=[('Prokaryote', 'Prokaryote'), ('Eukaryote', 'Eukaryote')], render_kw={"id": "cell_type"})
    submit = SubmitField("Submit")





@app.route('/generate', methods=['GET', 'POST'])
def generate():
    email = None
    gene = None
    trigger = None
    rna = None
    cell_type = None

    input_form = InputForm()
    if input_form.validate_on_submit():
        try:
            trigger = input_form.trigger.data.upper()
            gene = input_form.gene.data.upper()
        except:
            pass
        try:
            rna = input_form.rna.data.upper()
        except:
            pass

        # Get the data from the form
        email = input_form.email.data
        gene = input_form.gene.data
        trigger = input_form.trigger.data
        mrna = input_form.rna.data
        cell_type = input_form.cell_type.data


        s_email = str(email) if email else "EMPTY"
        s_gene = str(gene) if gene else "EMPTY"
        s_trigger = str(trigger) if trigger else "EMPTY"
        s_mrna = str(mrna) if mrna else "EMPTY"
        s_cell_type = str(cell_type) if cell_type else "EMPTY"


        # Run the scoring script in different process
        subprocess.Popen(["python3", "-m", "debugpy", "--listen", str(random.randint(10000, 50000)), "run_scoring.py",
                          s_gene, s_trigger, s_mrna, s_cell_type, s_email])

        # Clear the form
        input_form.gene.data = ''
        input_form.trigger.data = ''
        input_form.rna.data = ''
        input_form.email.data = ''
        input_form.cell_type.data = ''

    return render_template('generate.html', input_form=input_form, gene=gene, trigger=trigger, rna=rna, cell_type=cell_type)





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

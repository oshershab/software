import subprocess
import random
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
    input_form = InputForm()
    if input_form.validate_on_submit():
        # Get the data from the form and process it
        email = input_form.email.data
        gene = input_form.gene.data.upper()
        trigger = input_form.trigger.data.upper()
        mrna = input_form.rna.data.upper()
        cell_type = input_form.cell_type.data

        # Prepare the arguments for the subprocess
        args = [str(email), gene, trigger, mrna, cell_type]
        data = {'email': email, 'gene': gene, 'trigger': trigger, 'mrna': mrna, 'cell_type': cell_type}

        try:
            # Run the scoring script in a different process
            port = random.randint(10000, 50000)
            subprocess.Popen(["python3", "-m", "debugpy", "--listen", str(port), "run_scoring.py"] + args)
        except Exception as e:
            print(f"Failed to start subprocess: {e}")

        # Clear the form data
        input_form.email.data = ''
        input_form.gene.data = ''
        input_form.trigger.data = ''
        input_form.rna.data = ''
        input_form.cell_type.data = ''

        # Optional: add a success message or redirect
        return render_template('generate.html', input_form=input_form, success=True, data=data)

    return render_template('generate.html', input_form=input_form, success=False)





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

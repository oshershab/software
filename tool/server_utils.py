import os
import re
from flask import abort, current_app
from werkzeug.utils import secure_filename



# TODO: generic_process_file function, implement with biopython ? consider run time and memory usage.
def generic_process_file(uploaded_file):
    # create dictionary to store faste date from file
    fasta_dict = {}
    if uploaded_file:
        filename = secure_filename(uploaded_file.filename)
        file_ext = os.path.splitext(filename)[1]
        if file_ext not in current_app.config['UPLOAD_EXTENSIONS']:
            abort(400)
        else:
            header = ''
            for line in uploaded_file:
                try:
                    line = line.decode('utf-8').strip()
                    if line.startswith('>'):
                        header = line
                        fasta_dict[header] = ''
                        continue
                    else:
                        line = line.upper()
                        if validate_sequence(line):
                            fasta_dict[header] += line.strip().upper()
                except ValueError as e:
                    return str(e)

            if not fasta_dict:
                return None
    return fasta_dict



def validate_sequence(sequence):
    if not re.search(r"^[ACGT]", sequence):
        raise ValueError("Invalid sequence: must contain only A, C, G, or T nucleotides.")

    return sequence




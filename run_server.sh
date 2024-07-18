export FLASK_APP=ProTech.py
export FLASK_ENV=development

flask run --host 0.0.0.0 -p 80 > /tmp/server_log.txt 2> /tmp/server_error.txt &
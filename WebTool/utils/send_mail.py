import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

def send_mail(email_adress: str, pdf_path: str):
    mail_content = '''Hello,
    Attached to this mail is a detailed PDF with the results from your switch generator job.
    Thank You for choosing TrigGate'''
    #The mail addresses and password
    sender_address = '2022igemtau@gmail.com'
    sender_pass = 'hbetlnbeubtdhwex'
    receiver_address = email_adress

    #Setup the MIME
    message = MIMEMultipart()
    message['From'] = sender_address
    message['To'] = receiver_address
    message['Subject'] = "Results from TrigGate's switch generator"   #The subject line
    #The body and the attachments for the mail
    message.attach(MIMEText(mail_content, 'plain'))

    # open the file in bynary
    binary_pdf = open(pdf_path, 'rb')

    payload = MIMEBase('application', 'octate-stream', Name=pdf_path)
    payload.set_payload((binary_pdf).read())

    # enconding the binary into base64
    encoders.encode_base64(payload)

    # add header with pdf name
    payload.add_header('Content-Decomposition', 'attachment', filename=pdf_path)
    message.attach(payload)

    #Create SMTP session for sending the mail
    session = smtplib.SMTP('smtp.gmail.com', 587) #use gmail with port
    session.starttls() #enable security

    session.login(sender_address, sender_pass) #login with mail_id and password
    text = message.as_string()
    session.sendmail(sender_address, receiver_address, text)
    session.quit()
    return "Mail sent"

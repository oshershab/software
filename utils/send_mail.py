# Email configuration
import os
import smtplib
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

SMTP_SERVER = 'smtp.gmail.com'
SMTP_PORT = 587
EMAIL_ADDRESS = 'oshershab@gmail.com'
EMAIL_PASSWORD = 'uqme nvqb wgcf cgwl'


def send_email_with_attachment(to_email, subject, body, files):
    # Create the email message
    msg = MIMEMultipart()
    msg['From'] = EMAIL_ADDRESS
    msg['To'] = to_email
    msg['Subject'] = subject

    # Attach the email body
    msg.attach(MIMEText(body, 'plain'))

    # Attach files
    for file in files:
        if os.path.isfile(file):
            #print(f"Attaching file: {file}")  # Debugging line
            with open(file, "rb") as attachment:
                part = MIMEBase('application', 'octet-stream')
                part.set_payload(attachment.read())
                encoders.encode_base64(part)
                part.add_header('Content-Disposition', f"attachment; filename= {os.path.basename(file)}")
                msg.attach(part)
        else:
            pass
            #print(f"File not found: {file}")  # Debugging line
    return msg


def send_email(msg):
    # Connect to the server
    server = smtplib.SMTP(SMTP_SERVER, SMTP_PORT)
    server.starttls()  # Secure the connection
    server.login(EMAIL_ADDRESS, EMAIL_PASSWORD)  # Login to the email account

    # Send the email
    server.send_message(msg)
    server.quit()
    print("Email sent successfully!")


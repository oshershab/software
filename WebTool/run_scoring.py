from utils.send_mail import send_mail
from utils.pdf_generator import generate_pdf
import sys





def main(trigger, gene, mrna, cell_type, email):
    ###logic###






    ## send results to client ##
    pipeline_results = " "
    generate_pdf(pipeline_results, f"./", order_by='')
    send_mail(email, f"{pipeline_results.job_id}.pdf")





if __name__ == '__main__':
    # parameters from flask form
    trigger = sys.argv[1]
    gene = sys.argv[2]
    email = sys.argv[3]
    cell_type = sys.argv[4]
    mrna = sys.argv[5]

    # send to logic
    main(trigger, gene, mrna, cell_type, email)


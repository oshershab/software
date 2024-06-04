import sys

def main(trigger, gene, mrna, cell_type, email):
    print("Processing...")
    # logic here
    print(f"Trigger: {trigger}, Gene: {gene}, mRNA: {mrna}, Cell Type: {cell_type}, Email: {email}")







    ## send results ##
    #pipeline_results = " "
    #generate_pdf()
    #send_mail()


if __name__ == '__main__':

    if len(sys.argv) != 6:  # Expecting five arguments plus the script name
        print("Invalid number of arguments.")
        sys.exit(1)


    # parameters from flask form
    email = sys.argv[1]
    gene = sys.argv[2]
    trigger = sys.argv[3]
    mrna = sys.argv[4]
    cell_type = sys.argv[5]

    # send to logic
    main(trigger, gene, mrna, cell_type, email)




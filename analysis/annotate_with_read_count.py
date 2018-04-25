import pandas as pd
import re
from Bio import Entrez
import xml.etree.ElementTree as ET
import time
import argparse

parser = argparse.ArgumentParser(description='Annotate read count csv with spot count')

parser.add_argument('--input_csv', type=str, nargs=1, required=True,
                    help='The location of the read count csv file (input)')
parser.add_argument('--email', type=str, nargs=1, required=True,
                    help='Email for Entrez.')
parser.add_argument('--max_errors', type=int, nargs=1, required=True,
                    help='Number of times to try the API')

args = parser.parse_args()

def get_srr_accession(df):
    """
    Extract the SRR number from the file_location field of the dataframe.

    """
    pattern = re.compile(r"[SDE]RR[0-9]+")

    file = df['file_location']

    return re.search(pattern, file).group(0)

def get_spot_count(df, email, max_errors):
    """
    Use the Entrez API to get the spot (read) count for that SRR number.
    """

    srr_acc = df['srr_acc']

    # get SRR id
    Entrez.email = email

    error_count =0

    while error_count < max_errors:

        try:

            handle = Entrez.esearch(db="sra",term=srr_acc)

            record = Entrez.read(handle)

            handle.close()

            srr_id = record["IdList"][0]

            break

        except:

            print ('error occured collecting ID ', srr_acc, error_count)

            time.sleep(10)

            error_count = error_count +1


    # get SRR summary

    error_count =0

    while error_count < max_errors:

        try:

            handle = Entrez.esummary(db='sra', id=srr_id)

            record = Entrez.read(handle)

            handle.close()

            my_xml = record[0]['Runs']

            # Parse XML
            xml_object = ET.fromstringlist(["<root>", my_xml, "</root>"])

            # Get total_spots (reads)
            for child in xml_object:

                if child.attrib['acc'] == srr_acc:

                    print (srr_acc, child.attrib['total_spots'])

                    return int(child.attrib['total_spots'])

        except:

            print ('error occured collecting spot count', srr_acc, error_count)

            time.sleep(10)

            error_count = error_count +1

    return 1

# Load CSV
df = pd.read_csv(args.input_csv[0], index_col=0)

#Collect spot/read count data
df['srr_acc'] = df.apply(get_srr_accession, axis=1)
df['total_spots'] = df.apply(get_spot_count, axis=1, args=[args.email[0], args.max_errors[0]])
df['total_spots'] = df['total_spots'].astype(int)


#normalise the last three columns i.e. the exact hits covering the ROI
df['norm_alpha_read_covers_snps_count_exact'] = df['alpha_read_covers_snps_count_exact'] / (df['total_spots']/ 1000000)
df['norm_alpha_dup_read_covers_snps_count_exact'] = df['alpha_dup_read_covers_snps_count_exact'] / (df['total_spots']/ 1000000)
df['norm_beta_read_covers_snps_count_exact'] = df['beta_read_covers_snps_count_exact'] / (df['total_spots']/ 1000000)

output_file_name = '.'.join(args.input_csv[0].split('.')[:-1])+'.normal.csv'

df.to_csv(output_file_name)

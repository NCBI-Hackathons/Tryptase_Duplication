import argparse
import os



parser = argparse.ArgumentParser(description='Create a CSV file for further processing')

parser.add_argument('--bam_folder', type=str, nargs=1, required=True,
                    help='The directory containing the bam files.')
parser.add_argument('--source', type=str, nargs=1, required=True,
                    help='The source of the bam files e.g. AML')
parser.add_argument('--roi_start', type=int, nargs=1, required=True,
                    help='The 0 based start position of ROI')
parser.add_argument('--roi_end', type=int, nargs=1, required=True,
                    help='The 0 based end position of ROI')
parser.add_argument('--email', type=str, nargs=1, required=True,
                    help='Email for Entrez.')
parser.add_argument('--max_errors', type=int, nargs=1, required=True,
                    help='Number of times to try the API')
parser.add_argument('--file_comment', type=str, nargs=1, required=True,
                    help='text to add to file name')
args = parser.parse_args()
print (args)


#run create_read_count

output_name = '{source}.{file_comment}.csv'.format(
                    source=args.source[0],
                    file_comment=args.file_comment[0])


command = ("python create_read_count_df.py --bam_folder {bam_folder} --source {source} "
            "--roi_start {start} --roi_end {end} --output_name {output_name}").format(
                        bam_folder = args.bam_folder[0],
                        source = args.source[0],
                        start = args.roi_start[0],
                        end = args.roi_end[0],
                        output_name =output_name)


print (command)
os.system(command)


#annotate and normalise
command = ("python annotate_with_read_count.py --input_csv {input} "
            "--email {email} --max_errors {max_errors}").format(input = output_name,
            email=args.email[0],
            max_errors=args.max_errors[0])

print (command)
#os.system(command)

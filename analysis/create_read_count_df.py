import pysam
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Create a CSV file for further processing')

parser.add_argument('--bam_folder', type=str, nargs=1, required=True,
                    help='The directory containing the bam files.')
parser.add_argument('--source', type=str, nargs=1, required=True,
                    help='The source of the bam files e.g. AML')
parser.add_argument('--roi_start', type=int, nargs=1, required=True,
                    help='The 0 based start position of ROI')
parser.add_argument('--roi_end', type=int, nargs=1, required=True,
                    help='The 0 based end position of ROI')
parser.add_argument('--output_name', type=str, nargs=1, required=True,
                    help='The name of the output file')
args = parser.parse_args()

print (args)


def get_files(directory_path):
    """
    Return a list containing all the files matching a pattern in \
    the given directory

    """

    print ('{directory_path}/*sam.sorted.bam'.format(
                    directory_path=directory_path))
    file_list  = glob.glob('{directory_path}/*sam.sorted.bam'.format(
                    directory_path=directory_path))

    return file_list

def get_read_count(df):
    """
    Get the read count in a specific bam file.

    """

    return int(pysam.view("-c", df['file_location']))

def get_transcript_count(df, transcript_name):

    """
    Return the count of the transcript i.e how many of each of \
    the three transcripts are in there.
    Note - Does not look at the quality of the alignment.
    """

    sam_file_location = df['file_location']

    samfile = pysam.AlignmentFile(sam_file_location, "rb")

    count = 0

    for read in samfile:

        if read is not None and read.reference_name == transcript_name:

            count = count + 1

    return count

def get_zero_edit_distance_count(df, transcript_name):

    """
    For a goven transcript e.g. Alpha_GEX_64k_HEX get how many of \
    the matches have an NM tag of zero
    i.e. the match was exact.

    """

    sam_file_location = df['file_location']

    samfile = pysam.AlignmentFile(sam_file_location, "rb")

    count = 0

    for read in samfile:

        if read is not None and read.reference_name == transcript_name:

            edit_distance = read.get_tag('NM')

            if edit_distance == 0:

                count = count + 1

    return count

def get_transcript_read_count_filtered(df, transcript_name, start, end):

    """
    Count hits which cover the bit of the reference we are interested in.

    start = 0 based position in reference the alignment must start on or before.
    end = 0 based position that the alignment must end on or before.

    """

    sam_file_location = df['file_location']

    samfile = pysam.AlignmentFile(sam_file_location, "rb")

    iter = samfile.fetch(transcript_name, start, end)

    count =0

    for read in iter:

        if read.reference_start <=start and read.reference_end >= end:

            count = count +1

    return count

def get_transcript_read_count_filtered_exact(df, transcript_name, start, end):

    """
    Count hits which cover the bit of the reference we are interested in \
    and which are exact.

    That is - do they cross position 0 - 44 of the transcript \
    (given by transcipt_name) and \
    have an edit distance from the reference of 0.

    start = 0 based position in reference the alignment must start on or before.
    end = 0 based position that the alignment must end on or before.

    """

    sam_file_location = df['file_location']

    samfile = pysam.AlignmentFile(sam_file_location, "rb")

    iter = samfile.fetch(transcript_name, start, end)

    count =0

    for read in iter:

        if read.reference_start <=start and read.reference_end >=end:

            edit_distance = read.get_tag('NM')

            if edit_distance == 0:

                count = count + 1

    return count


results = get_files(args.bam_folder[0])
file_names = [x.split('/')[len(x.split('/'))-1] for x in results]

#Create Dataframe
df = pd.DataFrame(index=file_names)

df['file_location'] = results
df['source'] = args.source[0]


#Add total-count column
df['alignment_count'] = df.apply(get_read_count, axis=1)
print ('Added total alignment count')



#Count the number of reads aligned to each reference transcript
df['alpha_wt_count'] = df.apply(get_transcript_count, axis=1, args=['Alpha_GEX_64k_HEX'])
df['alpha_dup_count'] = df.apply(get_transcript_count, axis=1, args=['Alpha_GEX_79k_dup_FAM'])
df['beta_count'] = df.apply(get_transcript_count, axis=1, args=['BETA_new_GEX_FAM'])
print ('Added alignment count for each transcript')


#Number of exact read alignments for each reference transcript
df['alpha_wt_zero_edit_count'] = df.apply(get_zero_edit_distance_count, axis=1, args=['Alpha_GEX_64k_HEX'])
df['alpha_dup_zero_edit_count'] = df.apply(get_zero_edit_distance_count, axis=1, args=['Alpha_GEX_79k_dup_FAM'])
df['beta_zero_edit_count'] = df.apply(get_zero_edit_distance_count, axis=1, args=['BETA_new_GEX_FAM'])
print ('Added exact alignment count for each transcript')

# Number of alignments that span our area of interest

start = args.roi_start[0]
end = args.roi_end[0]

df['alpha_read_covers_snps_count'] = df.apply(get_transcript_read_count_filtered,
                                              axis=1,
                                              args=['Alpha_GEX_64k_HEX', start,end])
df['alpha_dup_read_covers_snps_count'] = df.apply(get_transcript_read_count_filtered,
                                              axis=1,
                                              args=['Alpha_GEX_79k_dup_FAM', start,end])
df['beta_read_covers_snps_count'] = df.apply(get_transcript_read_count_filtered,
                                              axis=1,
                                              args=['BETA_new_GEX_FAM', start,end])
print ('Added alignment count for each transcript within ROI')


# Number of alignments that span our area of interest and are exact
df['alpha_read_covers_snps_count_exact'] = df.apply(get_transcript_read_count_filtered_exact,
                                              axis=1,
                                              args=['Alpha_GEX_64k_HEX', start,end])

df['alpha_dup_read_covers_snps_count_exact'] = df.apply(get_transcript_read_count_filtered_exact,
                                              axis=1,
                                              args=['Alpha_GEX_79k_dup_FAM', start,end])

df['beta_read_covers_snps_count_exact'] = df.apply(get_transcript_read_count_filtered_exact,
                                              axis=1,
                                              args=['BETA_new_GEX_FAM', start,end])

print ('Added exact alignment count for each transcript within ROI')

df.to_csv(args.output_name[0])
print ('Created CSV')

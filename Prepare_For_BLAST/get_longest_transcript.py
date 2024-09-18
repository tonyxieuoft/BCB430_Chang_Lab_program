import os
from Basic_Tools.lists_and_files import file_to_list


def get_longest_transcript(directory):
    """
    Return the f_read_name of the longest transcript within a specified directory.

    :param directory: a directory containing transcripts in fasta format, where
    each transcript is separated into exon chunks and the length of each
    transcript is in the fasta heading.
    :return:
    """
    # TODO how to replace this for blast?
    longest_transcript_length = 0
    longest_transcript = ""
    for transcript in os.listdir(directory):
        transcript_arr = file_to_list(os.path.join(directory, transcript))

        exon_section = 0
        sample_header = transcript_arr[0].split(" ")
        while len(sample_header[exon_section]) < 6 or \
                sample_header[exon_section][:6] != "genome":
            exon_section += 1
        exon_section += 1

        transcript_length = 0
        line_no = 0
        while line_no < len(transcript_arr):
            if transcript_arr[line_no][0] == ">":
                exon_endpoints = transcript_arr[line_no].split(" ")[exon_section].split("-")
                transcript_length += int(exon_endpoints[1]) - int(exon_endpoints[0]) + 1
            line_no += 1

        if transcript_length > longest_transcript_length:
            longest_transcript = transcript
            longest_transcript_length = transcript_length

    return longest_transcript


if __name__ == "__main__":

    directory = r"C:\Users\tonyx\Downloads\Balaenoptera musculus"
    #print(get_longest_transcript(directory))




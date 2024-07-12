from Basic_Tools.lists_and_files import file_to_list


def concatenate_exons(fasta_path: str) -> str:
    """
    Concatenates the exons pointed to in fasta_path into one transcript

    :param fasta_path: path containing exons
    :return: concatenated exons into one transcript, as a fasta string
    """

    # get the exons file into a list
    fasta_list = file_to_list(fasta_path)
    # this is the fasta title, split into sections
    title_sections = fasta_list[0].split(" ")

    # find the index of the genome section
    section_no = 0
    while len(title_sections[section_no]) < 6 or \
            title_sections[section_no][:6] != "genome":
        section_no += 1

    # the exon section is right after the genome section
    exon_section = section_no + 1

    # earliest index we have of the gene
    global_beginning = title_sections[exon_section].split("-")[0]

    sequence = ""
    last_range = None
    line_no = 0
    while line_no < len(fasta_list):

        # current exon boundary
        coverage = fasta_list[line_no].split(" ")[exon_section].split("-")

        line_no += 1
        # I need to know why this works lol
        if line_no < len(fasta_list):

            if fasta_list[line_no][0] == ">":
                # executes only when there's no fasta sequence after a heading, for
                # some reason -> the earlier one is if it's the last, and the latter
                # is if we see two fasta headers in a row
                missing_space = (int(coverage[1]) - int(coverage[0]) + 1)
                # missing space is just whatever is in the header
                sequence += "-"*missing_space
            elif last_range is None:
                # the first fasta sequence that's encountered
                # if there is a gap in the beginning, use '-'
                if global_beginning != 1:
                    sequence += (int(global_beginning) - 1)*"-"
                sequence += fasta_list[line_no]
                line_no += 1
            else:
                # track gaps between the current sequence and the last, fills with
                # '-' if that's the case
                in_between = int(coverage[0]) - int(last_range[1]) - 1
                filler = "-"*in_between
                sequence += filler + fasta_list[line_no]
                line_no += 1

        last_range = coverage

    # update fasta heading with beginning and end parameters of whole gene
    fasta_heading = ""
    for i in range(exon_section):
        if i == 0:
            fasta_heading = title_sections[i]
        else:
            fasta_heading += " " + title_sections[i]

    fasta_heading += " " + global_beginning + "-" + last_range[1] + "\n"

    return fasta_heading + sequence + "\n"


if __name__ == "__main__":
    fasta_path = r"C:\Users\tonyx\Downloads\blast_test_results1234567891011121314151617\rdh8\lamniformes\Carcharodon carcharias\Reference_XM_041177394.1.fas"
    concatenate_exons(fasta_path)

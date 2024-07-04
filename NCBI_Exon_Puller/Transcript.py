
class Transcript:

    def __init__(self, accession: str):

        self.accession = accession

        self.transcript_sequence = []
        self.cds = []

        self.cds_start = -1
        self.cds_end = 0

    def check_valid_cds(self, exon, cds) -> bool:

        if min(int(exon[0]), int(exon[1])) > min(int(cds[0]), int(cds[1])):
            return False

        if max(int(exon[0]), int(exon[1])) < max(int(cds[0]), int(cds[1])):
            return False

        return True





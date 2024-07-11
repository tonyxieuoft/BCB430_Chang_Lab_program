
class NoGenomicBeginningException(Exception):
    pass

class Transcript:

    def __init__(self, accession: str):

        self.accession = accession

        self.transcript_sequence = []
        self.cds = []

        self.cds_start = -1
        self.cds_end = 0

        self.exons_no_cds_space = 0

    def check_valid_cds(self, exon, cds) -> bool:

        if min(int(exon[0]), int(exon[1])) > min(int(cds[0]), int(cds[1])):
            return False

        if max(int(exon[0]), int(exon[1])) < max(int(cds[0]), int(cds[1])):
            return False

        return True

    def add_to_exons_no_cds_space(self, exon):
        self.exons_no_cds_space += abs(int(exon[1]) - int(exon[0])) + 1

    def add_cds(self, exon, cds):

        cds_range = abs(int(cds[1]) - int(cds[0]))
        if self.cds_start == -1:
            self.cds_start = abs(int(cds[0]) - int(exon[0])) + 1 + self.exons_no_cds_space
            self.cds.append((self.cds_start, self.cds_start + cds_range))
            self.cds_end = self.cds_start + cds_range
        else:
            self.cds.append((self.cds_end + 1, self.cds_end + 1 + cds_range))
            self.cds_end = self.cds_end + 1 + cds_range






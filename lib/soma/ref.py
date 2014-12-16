import sys

class Reference(object):

    # FASTA Specifications
    id = ""
    sequence = ""

    def __init__(self, filename):
        self.filename = filename
        self.load_object()

    def load_object(self):
        with open(self.filename, "r") as f:
            for g in f:
                if g.startswith(">"):
                    self.id = g.strip("\r\n").lstrip(">").strip()
                else:
                    self.sequence += g.strip("\r\n").strip().upper()

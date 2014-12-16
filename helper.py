import sys
import time
import re


class Helper(object):
    
    def get_time(self):
        return time.strftime("%Y-%m-%d %H:%M:%S")

    def get_other_base(self, base, indel=False):
        if indel:
            if base == "A":
                return ["T", "G", "C", "I", "D"]
            elif base == "T":
                return ["A", "G", "C", "I", "D"]
            elif base == "G":
                return ["A", "T", "C", "I", "D"]
            elif base == "C":
                return ["A", "T", "G", "I", "D"]
            elif base == "I":
                return ["A", "T", "G", "C", "D"]
            elif base == "D":
                return ["A", "T", "G", "C", "I"]
            else:
                return []
        else:
            if base == "A":
                return ["T", "G", "C"]
            elif base == "T":
                return ["A", "G", "C"]
            elif base == "G":
                return ["A", "T", "C"]
            elif base == "C":
                return ["A", "T", "G"]
            else:
                return []
    
    def convert_to_str(self, pc_genome_gt):
        return "%.2f" % (pc_genome_gt)

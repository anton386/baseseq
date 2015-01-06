import sys
import math
import re
import pprint as pp


import pysam

from helper import Helper


class Consensus(Helper):

    def __init__(self, bam, ref):
        self.bam = bam
        self.ref = ref
        self.consensus = {}
        self.consensus_genomes = {}     # 1-based
        self.consensus_posteriors = {}  # 1-based
        self.consensus_indels = {}      # 1-based
        self.consensus_coverages = {}   # 1-based
        self.reads_in_barcode = {}

        # Summary statistics
        self.no_of_reads_analyzed = 0
        self.no_of_complete_genomes = 0
        self.no_of_high_quality_genomes = 0

        # Constants
        self.log_table = self.generate_log_table()
        self.error_table = self.generate_error_table()
        self.base_error_table = dict()
        self.indel_error_table = dict()
        
        self.regex_cigar_num_base_to_skip = re.compile("[0-9]+")
        self.regex_cigar_alpha_operator = re.compile("[A-Z]")

    def generate_log_table(self):
        log_table = {}
        for i in range(1, 42):
            log_table[i] = -math.log10(1.00 - self.phred_to_pval(i))

        return log_table

    def generate_error_table(self):
        error_table = {}
        for i in range(1, 42):
            error_table[i] = -math.log10((self.phred_to_pval(i) / 5))

        return error_table

    def generate_base_error_table(self, no_of_indels, e_indel=45):
        p_indel = self.phred_to_pval(e_indel)
        l_indel = self.pval_to_log(p_indel)
        error_table = {}
        for i in range(1, 42):
            p_base = ((self.phred_to_pval(i) - no_of_indels*p_indel)/3)  # uniform base error for indel
            l_base = self.pval_to_log(p_base)
            error_table[i] = {"A": l_base, "T": l_base, "G": l_base, "C": l_base,
                              "I": l_indel, "D": l_indel}

        return error_table

    def generate_indel_error_table(self, no_of_indels, e_indel=45):
        p_indel = self.phred_to_pval(e_indel)
        l_indel = self.pval_to_log(p_indel)
        error_table = {}
        for i in range(1, 42):
            p_base = ((self.phred_to_pval(i) - ((no_of_indels*p_indel)-1))/4)  # uniform base error for indel
            l_base = self.pval_to_log(p_base)
            error_table[i] = {"A": l_base, "T": l_base, "G": l_base, "C": l_base,
                              "I": l_indel, "D": l_indel}

        return error_table

    def phred_to_pval(self, phred):
        return 10**(-phred/10.0)

    def pval_to_phred(self, pval):
        return -10.0*(math.log10(pval))

    def log_to_pval(self, log):
        return 10**(-log)

    def pval_to_log(self, pval):
        return -math.log10(pval)

    
    def build(self, phred=33, q_threshold=0, q_indel=40,
              high_quality_threshold=0.15):
        #TODO use EM to estimate haplotype frequencies
        #TODO haplotype imputation

        bam = pysam.Samfile(self.bam, "rb")

        current_barcode = ""
        consensus_matrix = self.create_log_matrix()
        consensus_quality = self.create_quality_matrix()
        consensus_coverage = self.create_coverage_matrix()
        no_of_barcodes = 0
        
        for rindex, samy in enumerate(bam):

            self.check_status_of_reads(rindex)
            
            barcode = samy.qname.split("|")[0].replace("BC:", "")
            
            if (current_barcode != "") and (current_barcode != barcode):
                #DEBUG print self.reads_in_barcode[current_barcode]
                #DEBUG print current_barcode
                #DEBUG pp.pprint(consensus_matrix)
                
                # create consensus_genome
                # calculate consensus_posterior
                # tally coverage
                
                self.calculate_score(consensus_matrix,
                                     consensus_quality)
                
                self.create_consensus(consensus_matrix,
                                      consensus_coverage, current_barcode,
                                      high_quality_threshold=high_quality_threshold)
                self.tally_coverage(consensus_coverage, current_barcode)

                self.check_status_of_barcodes(no_of_barcodes)
                no_of_barcodes += 1
                
                current_barcode = barcode
                consensus_matrix = self.create_log_matrix()
                consensus_quality = self.create_quality_matrix()
                consensus_coverage = self.create_coverage_matrix()
                
            
            current_barcode = barcode

            try:
                self.reads_in_barcode[current_barcode] += 1
            except KeyError:
                self.reads_in_barcode[current_barcode] = 1
            
            numeric = map(int, self.regex_cigar_num_base_to_skip.findall(samy.cigarstring))
            operator = self.regex_cigar_alpha_operator.findall(samy.cigarstring)

            # BETA adjust CIGAR
            if "I" in operator:
                index_of_I = []
                for n, op in enumerate(operator):
                    if op == "I":
                        index_of_I.append(n)
                for n in index_of_I:
                    if operator[n-1] == "M":
                        numeric[n-1] = numeric[n-1] - 1  # subtract 1 base from match
                        numeric[n] =  numeric[n] + 1  # add 1 base to insertion
                adjusted = True
            else:
                adjusted = False
            
            start = 0
            ref_start = 0
            for num, (c, op) in enumerate(zip(numeric, operator)):
                if op == "I":
                    position = int(samy.pos+1) + ref_start
                    #adj_pos = start - 1  # adjusted to include ref, not because pysam is 0-based
                    indel = (samy.seq[start:start+c], c, op)
                    allele = samy.seq[start:start+c]
                    
                    if "N" in allele:  # ignore allele if there is N in the insertion
                        continue
                    
                    q = q_indel
                    try:
                        consensus_matrix[position]
                        consensus_quality[position]
                        consensus_coverage[position]

                        try:
                            consensus_matrix[position]["INS"]

                            # remove INS and replace with new allele
                            del consensus_matrix[position]["INS"]
                            consensus_matrix[position][allele] = 0
                            
                        except KeyError:
                            consensus_matrix[position][allele] = 0

                        try:
                            consensus_quality[position]["INS"]

                            # remove INS and replace with new allele
                            del consensus_quality[position]["INS"]
                        except KeyError:
                            try:
                                consensus_quality[position][allele].append(q)
                            except KeyError:
                                consensus_quality[position][allele] = [q]
                        
                        consensus_coverage[position] += 1
                    
                    except KeyError:
                        # because of insertions after the genome
                        pass
                    
                    if adjusted and (num in index_of_I):
                        start += c
                        ref_start += 1  # move one ahead
                    else:
                        start += c
                        
                    
                elif op == "D":
                    position = int(samy.pos+1) + ref_start
                    #adj_pos = start - 1
                    indel = ("-" * c, c, op)

                    q = q_indel
                    try:
                        consensus_matrix[position]
                        consensus_quality[position]
                        consensus_coverage[position]

                        for no in range(c):
                            consensus_quality[position + no]["D"].append(q)
                            consensus_coverage[position + no] += 1
                        
                    except KeyError:
                        # because of deletions after the genome
                        pass
                    
                    ref_start += c
                elif op == "S":
                    start += c
                elif op == "H":
                    pass
                elif op == "M":

                    #if adjusted and (num+1 in index_of_I):
                    #    print c
                    #    print samy.pos + 1 + ref_start
                    #    print samy.seq[start:start+c]
                        
                    for no, (base, quality) in enumerate(zip(samy.seq[start:start+c], samy.qual[start:start+c])):
            
                        q = ord(quality) - phred
                
                        if q >= q_threshold:
                            try:

                                position = samy.pos + 1 + ref_start + no
                            
                                # add quality to array
                                consensus_quality[position][base].append(q)

                                # increment in consensus
                                consensus_coverage[position] += 1
                            except KeyError:
                                # because of N
                                pass
                    
                    ref_start += c
                    start += c
        
        # final step
        self.create_consensus(consensus_matrix, consensus_coverage,
                              current_barcode, high_quality_threshold=high_quality_threshold)
        self.tally_coverage(consensus_coverage, current_barcode)

        
        # output statistics
        sys.stderr.write("[%s] Number of Reads Analyzed: %s\n" % (self.get_time(),
                                                                  sum(self.reads_in_barcode.values())))
        sys.stderr.write("[%s] Number of BarCodes: %s\n" % (self.get_time(),
                                                            len(self.reads_in_barcode.keys())))
        sys.stderr.write("[%s] Number of High Quality BarCodes: %s\n" % (self.get_time(),
                                                                         len(self.consensus_genomes)))
        sys.stderr.write("[%s] Number of High Quality Genomes: %s\n" % (self.get_time(),
                                                                        self.no_of_high_quality_genomes))
        sys.stderr.write("[%s] Number of Complete Genomes: %s\n" % (self.get_time(),
                                                                    self.no_of_complete_genomes))
    

    def check_status_of_reads(self, rindex, interval=1000000):
        if rindex % interval == 0:
            sys.stderr.write("[%s] Currently at Read No: %s\n" % (self.get_time(),
                                                                  rindex))

    
    def check_status_of_barcodes(self, no_of_barcodes, interval=1000):
        if no_of_barcodes % interval == 0:
            sys.stderr.write("[%s] Currently at Barcode No: %s\n" % (self.get_time(),
                                                                     no_of_barcodes))

    
    def calculate_score(self, consensus_matrix, consensus_quality):
        for position, var in consensus_quality.items():

            no_of_indels = 0
            for allele in var.keys():
                if len(allele) > 1:
                    no_of_indels += 1
                else:
                    if allele == "D":
                        no_of_indels += 1

            # initialize base and indel error table
            self.base_error_table = self.generate_base_error_table(no_of_indels)
            self.indel_error_table = self.generate_indel_error_table(no_of_indels)
            all_other_base = self.get_all_other_base(var.keys())
            
            for allele, qual in var.items():
                if len(allele) == 1:  # either base or deletion
                    if allele == "D":
                        for q in qual:
                            consensus_matrix[position][allele] += self.log_table[q]

                            for other in all_other_base[allele]:
                                if len(other) > 1:
                                    consensus_matrix[position][other] += self.indel_error_table[q]["I"]
                                else:
                                    consensus_matrix[position][other] += self.indel_error_table[q][other]
                    else:
                        for q in qual:
                            consensus_matrix[position][allele] += self.log_table[q]

                            for other in all_other_base[allele]:
                                if len(other) > 1:
                                    consensus_matrix[position][other] += self.base_error_table[q]["I"]
                                else:
                                    consensus_matrix[position][other] += self.base_error_table[q][other]
                else:  # insertion
                    for q in qual:
                        consensus_matrix[position][allele] += self.log_table[q]

                        for other in all_other_base[allele]:
                            consensus_matrix[position][other] += self.indel_error_table[q][other]


    def create_consensus(self, consensus_matrix, coverage, barcode,
                         high_quality_threshold=0.15, minimum_reads_per_site=4):
        #if barcode == "AAAAAAGAACAGCCTGCATC":
        #pp.pprint(consensus_matrix)
        
        consensus_genome = {}
        consensus_posterior = {}
        
        for pos, matrix in sorted(consensus_matrix.items(), key=lambda q: q[0]):
            min_value = 10000000.0
            min_base = "N"
            for base, qual in matrix.items():
                if qual < min_value:
                    min_value = qual
                    min_base = base

            if min_value == 0.00:
                min_base = "N"
                min_posterior = 0
            else:
                min_posterior = self.calculate_max_posterior(min_value, matrix.values())

            if coverage[pos] < minimum_reads_per_site:
                min_base = "N"
                min_posterior = 0
            
            consensus_genome[pos] = min_base
            consensus_posterior[pos] = min_posterior

        
        #DEBUG print "".join(consensus_genome)
        #DEBUG print " ".join(map(str, consensus_posterior))
        no_of_missing_base = consensus_genome.values().count("N")
        percent_missing_base = float(no_of_missing_base) / float(len(self.ref.sequence))
        if no_of_missing_base == 0:
            self.no_of_complete_genomes += 1
        if percent_missing_base < high_quality_threshold:
            self.no_of_high_quality_genomes += 1

        
        # assign and store to self by barcode
        self.consensus_genomes[barcode] = consensus_genome
        self.consensus_posteriors[barcode] = consensus_posterior
        
    
    def infer_consensus(self, consensus, consensus_id="chrConsensus"):
        
        f_out = open(consensus, "w")
        
        output = []
        
        # consensus - count method

        no_of_barcodes = len(self.consensus_genomes)
        length_of_genome = len(self.ref.sequence)
        
        genomes = self.consensus_genomes.values()  # 0-based

        # work on indels first
        # accumulate insertions for simple consensus counting
        consensus_indels = {}  # 1-based
        self.inferred_consensus_indels = {}
        for barcode, v in self.consensus_indels.items():  # 1-based
            for pos, cindel in v.items():
                try:
                    consensus_indels[pos]
                except KeyError:
                    consensus_indels[pos] = {}
                
                try:
                    consensus_indels[pos][cindel] += 1
                except KeyError:
                    consensus_indels[pos][cindel] = 1

        # next, work on the rest
        for i in range(length_of_genome):

            consensus = {"A": 0, "T": 0, "G": 0, "C": 0, "I": 0, "D": 0, "N": 0}
            for j in range(no_of_barcodes):
                consensus[genomes[j][i]] += 1

            consensus_base = "N"
            consensus_count = 0
            for b, c in consensus.items():
                if c > consensus_count:
                    consensus_base = b
                    consensus_count = c

            # get the correct indel, otherwise, change consensus to "N"
            # since there is no indel here
            if consensus_base == "I":
                max_indel = None
                max_count = 0
                try:
                    for k, l in consensus_indels[i+1].items():  # i is 0-based
                        if max_count < l:
                            max_indel = k
                            max_count = l
                    self.inferred_consensus_indels[i+1] = k
                except KeyError:
                    consensus_base = "N"
                
            
            output.append(consensus_base)

        
        self.inferred_consensus = "".join(output)
        self.inferred_fconsensus = self.adjust_consensus(self.inferred_consensus_indels)
        
        f_out.write(">%s\n%s\n" % (consensus_id, self.inferred_consensus))

        f_out.close()


    def adjust_consensus(self, indels):
        consensus = {}
        for pos, base in enumerate(self.inferred_consensus):
            if base == "D":
                pass
            elif base == "I":
                consensus[pos+1] = indels[pos+1]  # pos is 0-based
            else:
                consensus[pos+1] = base
        
        ordered = []
        for k, v in sorted(consensus.items(), key=lambda q: q[0]):
            ordered.append(v)
        
        return "".join(ordered)
        

    def output_consensus_genomes(self, out):
        f_out = open(out, "w")

        for barcode, genome in self.consensus_genomes.items():
            f_out.write(">%s\n%s\n" % (barcode, "".join(genome)))

        f_out.close()

    def output_haplotype_distribution(self, out):
        f_out = open(out, "w")
        
        total = 0
        freq_distribution = {}
        for genome in self.consensus_genomes.values():
            haplotype = self.get_haplotype(genome)
            try:
                freq_distribution[haplotype] += 1
            except KeyError:
                freq_distribution[haplotype] = 1

            total += 1

        for haplotype, counts in sorted(freq_distribution.items(), key=lambda q: q[1], reverse=True):
            f_out.write("%s\t%s\t%s\n" % (haplotype, counts, float(counts)/float(total)))

        f_out.close()

    def output_vcf(self, out):
        f_out = open(out, "w")
        
        variants = {}
        total = 0
        for barcode, cg in self.consensus_genomes.items():
            total += 1
            for no, (x, y) in enumerate(zip(self.inferred_consensus, cg)):
                if x != y:
                    if x == "N":
                        continue
                    elif y != "N":
                        try:
                            variants[no+1]
                        except KeyError:
                            variants[no+1] = {}

                        try:
                            variants[no+1][y] += 1
                        except KeyError:
                            variants[no+1][y] = 1

        f_out.write("##fileformat=VCFv4.1")
        f_out.write("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("CHROM", "POS",
                                                           "ID", "REF",
                                                           "ALT", "QUAL",
                                                           "FILTER", "INFO"))
        for k, v in sorted(variants.items()):
            for base, freq in sorted(v.items()):
                f_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("consensus", k,
                                                                  ".", self.inferred_consensus[k-1],
                                                                  base, ".",
                                                                  ".",
                                                                  ";".join(["AF=%s" % (float(freq)/float(total)),
                                                                            "AN=%s" % (freq,)])))
        
        f_out.close()

    def output_consensus_coverage(self, f_out):

        f_out.write("##ConsensusCoverage\n")
        f_out.write("##" + "\t".join(["barcode",
                                      "mean",
                                      "mean_gt_zero",
                                      "median",
                                      "median_gt_zero",
                                      "pc_genome_gt"]) + "\n")
        for barcode, values in sorted(self.consensus_coverages.items(), key=lambda q: q[0]):
            f_out.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n" % (barcode,
                                                              values["mean"],
                                                              values["mean_gt_zero"],
                                                              values["median"],
                                                              values["median_gt_zero"],
                                                              " | ".join(map(self.convert_to_str,
                                                                            values["pc_genome_gt"]))))
    
    
    def get_haplotype(self, genome):
        haplotype = []
        for no, (i, j) in enumerate(zip(self.inferred_consensus, genome)):
            if i != j:
                if i == "N":
                    continue
                elif j != "N":
                    haplotype.append((no+1, j))
        
        return tuple(haplotype)

        
    def calculate_max_posterior(self, min_value, matrix):
        total_pval = 0.00
        for p in matrix:
            total_pval += self.log_to_pval(p)

        try:
            total_log = self.pval_to_log(total_pval)
            
            return int(10 * ((sum(matrix) - min_value) - total_log))
        
        except ValueError:
            # overflow error
            return 999999
    
    def create_log_matrix(self):
        matrix = {}
        for i in range(len(self.ref.sequence)):
            matrix[i+1] = {"A": 0, "T": 0, "G": 0, "C": 0, "INS": 0, "D": 0}

        return matrix

    def create_quality_matrix(self):
        matrix = {}
        for i in range(len(self.ref.sequence)):
            matrix[i+1] = {"A": [], "T": [], "G": [], "C": [], "INS": [], "D": []}

        return matrix

    def create_coverage_matrix(self):
        matrix = {}
        for i in range(len(self.ref.sequence)):
            matrix[i+1] = 0

        return matrix

    def tally_coverage(self, consensus_coverage, barcode):
        length_of_ref = len(consensus_coverage)
        cov = consensus_coverage.values()

        # mean
        mean = float(sum(cov)) / length_of_ref

        # median
        if (length_of_ref % 2) == 0:
            median = (cov[(length_of_ref/2)-1] + cov[(length_of_ref/2)+1]) / 2
        else:
            median = cov[(length_of_ref/2)]

        # mean gt 0
        sum_gt_zero = 0
        total_gt_zero = 0
        list_gt_zero = []
        for c in cov:
            if c > 0:
                sum_gt_zero += c
                total_gt_zero += 1
                list_gt_zero.append(c)
        mean_gt_zero = float(sum_gt_zero)/float(total_gt_zero)

        # median gt 0
        tgtzm = len(list_gt_zero)
        if (tgtzm % 2) == 0:
            median_gt_zero = (list_gt_zero[(tgtzm/2)-1] + list_gt_zero[(tgtzm/2)+1]) / 2
        else:
            median_gt_zero = list_gt_zero[(tgtzm/2)]

        # percent of genome greater than X
        pc_genome_gt = []
        for i in range(0, 31):
            total = 0
            for c in cov:
                if c > i:
                    total += 1
            pc_genome_gt.append(float(total)/length_of_ref)
        
        self.consensus_coverages[barcode] = {"median": median,
                                             "mean": mean,
                                             "median_gt_zero": median_gt_zero,
                                             "mean_gt_zero": mean_gt_zero,
                                             "pc_genome_gt": pc_genome_gt}

import sys
import re

import matplotlib
matplotlib.use("Agg")
import pysam
import matplotlib.pyplot as pyplot
import networkx as nx

from helper import Helper


class BarCode(Helper):

    
    def __init__(self, bam):
        self.bam = pysam.Samfile(bam, "rb")
        self.barcodes = []
        self.barcode_to_read = {}
        self.read_to_barcode = {}
        
        self.barcode_graph = None
        self.regex_cigar_num_base_to_skip = re.compile("[0-9]+")
        self.regex_cigar_alpha_operator = re.compile("[A-Z]")

    
    def cluster_barcodes(self, k=1):
        total_barcodes = 0
        self.barcode_graph = nx.Graph()
        
        for key in self.barcode_to_read.keys():
            for no, base in enumerate(key):
                for newbase in self.get_other_base(base):
                    newbarcode = list(key)
                    newbarcode[no] = newbase
                    newbarcode = "".join(newbarcode)
                    try:
                        self.barcode_to_read[newbarcode]
                        self.barcode_graph.add_edge(key, newbarcode)
                    except KeyError:
                        pass

            total_barcodes += 1

            if total_barcodes % 100000 == 0:
                sys.stderr.write("[%s] Total BarCodes: %s\n" % (self.get_time(), total_barcodes))

        
        clusters = nx.components.connected_components(self.barcode_graph)
        
        clusters_data = []

        total_barcodes_in_clusters = 0
        for no, cluster in enumerate(clusters):
            total_reads_in_barcodes = 0
            for bc in cluster:
                total_reads_in_barcodes += len(self.barcode_to_read[bc])
            print no, len(cluster), total_reads_in_barcodes
            total_barcodes_in_clusters += len(cluster)
            clusters_data.append(len(cluster))

        total_clusters = len(clusters)
        total_singletons = total_barcodes - total_barcodes_in_clusters
        print "Total Clusters: %s" % (total_clusters, )
        print "Total Singletons: %s" % (total_singletons, )

        pyplot.hist(clusters_data, bins=500)
        pyplot.savefig("clusters.png", dpi=200)

    
    def filter_and_export(self, list_of_ids, out, paired_end=False, export="fastq"):
        # output using paired_end mode ensures that the sequences are not treated as single end

        if export == "fastq":
            out = open(self.out, "w")
        elif export == "bam":
            out = pysam.Samfile(out, "wb", template=self.bam)
        
        new_list = []
        for ids in list_of_ids:
            new_list.append(int(ids.split(".")[1]))
        
        iter_of_ids = iter(sorted(new_list))
        read_id = iter_of_ids.next()

        total = 0
        total_reads = 0
        for samy in self.bam:

            if total_reads % 1000000 == 0:
                sys.stderr.write("[%s] Total Reads iterated: %s\n" % (self.get_time(), total_reads))

            bam_read_id = int(samy.qname.split(".")[1])
            
            if bam_read_id < read_id:
                total_reads += 1
                continue
            elif bam_read_id > read_id:
                try:
                    read_id = iter_of_ids.next()
                except StopIteration:
                    break
            else:

                if export == "fastq":
                    if paired_end:
                        pass
                    else:
                        out.write("\n".join(["@S7_%s" % (total, ),
                                             samy.query,
                                             "+",
                                             samy.qqual]) + "\n")
                elif export == "bam":
                    out.write(samy)

                total += 1

                total_reads += 1
                

    def load_barcodes(self, f_barcodes, reads_in_barcode_threshold=300):
        self.barcodes = []
        self.barcode_to_read = {}
        self.read_to_barcode = {}
        
        no = 0
        with open(f_barcodes, "r") as f:
            for line in f:
                data = line.strip("\r\n").split("\t")
                reads = data[1].split(",")
                
                if len(reads) >= reads_in_barcode_threshold:

                    self.barcodes.append(data[0])
                    
                    for read_id in data[1].split(","):
                        self.read_to_barcode[read_id] = no

                    no += 1

                    # just for now
                    self.barcode_to_read[data[0]] = data[1].split(",")


    def output_reads_in_barcode_distribution(self, f_out):
        f_out.write("##ReadsInBarcodeDistribution\n")
        f_out.write("##no_of_reads\tcounts\n")
        
        read_distribution = {}
        for barcode, reads in self.barcode_to_read.items():
            try:
                read_distribution[len(reads)] += 1
            except KeyError:
                read_distribution[len(reads)] = 1

        #TODO histogram with bins
        #for i in range(min(read_distribution.keys()), max(read_distribution.keys())+1):
        for i in sorted(read_distribution.keys()):
            try:
                f_out.write("%s\t%s\n" % (i, read_distribution[i]))
            except KeyError:
                f_out.write("%s\t%s\n" % (i, 0))

    
    def simple_approach(self):
        bcid = 0
        total = 0
        bcplaceholder = {}
        for samy in self.bam:
            if samy.is_proper_pair and samy.is_read2:
                # read is in a proper pair and is read2 in sequencing
                # barcodes are attached right next to P7

                numeric = map(int, self.regex_cigar_num_base_to_skip.findall(samy.cigarstring))
                operator = self.regex_cigar_alpha_operator.findall(samy.cigarstring)
                matched = False
                
                for c, op in zip(numeric, operator):
                    if op == "M":
                        matched = True
                    if op == "S" and matched:
                        if c >= 20:
                            #uniA = samy.seq[-c:-20][-1]
                            barcode = samy.seq[-20:]
                            try:
                                bcindex = bcplaceholder[barcode]
                                self.barcode_to_read[barcode].append(samy.qname)
                                
                            except KeyError:
                                self.barcodes.append(barcode)
                                self.barcode_to_read[barcode] = [samy.qname]
                                
                                
                                bcplaceholder[barcode] = bcid
                                bcindex = bcid

                                # increment
                                bcid += 1

                            self.read_to_barcode[samy.qname] = bcindex

                            total += 1

                            if total % 1000000 == 0:
                                sys.stdout.write("[%s] Loading Barcodes: %s\n" % (self.get_time(), total))
                            #sys.stdout.write(uniA + "\n")
                        else:
                            barcode = samy.seq[-c:]

    
    def sort_and_rewrite_bam(self, f_out):
        
        out = pysam.Samfile(f_out, "wb", template=self.bam)

        for samy in self.bam:
            try:
                barcode = self.barcodes[self.read_to_barcode[samy.qname]]
                samy.qname = "BC:%s|" % (barcode,) + samy.qname
                out.write(samy)
            except KeyError:
                continue
        
        out.close()

    def split_bam_into_barcodes(self, prefix=""):

        if not prefix:
            prefix = ""

        current_barcode = ""
        
        for samy in self.bam:

            barcode = samy.qname.split("|")[0].replace("BC:", "")

            if current_barcode == "":
                
                current_barcode = barcode
                f_out = self.get_split_bam_name(prefix, current_barcode)
                out = pysam.Samfile(f_out, "wb", template=self.bam)
                
            elif (current_barcode != "") and (current_barcode != barcode):

                out.close()

                current_barcode = barcode
                f_out = self.get_split_bam_name(prefix, current_barcode)
                out = pysam.Samfile(f_out, "wb", template=self.bam)
    
            out.write(samy)

        out.close()

    def get_split_bam_name(self, prefix, barcode):
        
        return prefix + barcode + ".bam"

    
    def write_barcodes(self, out):
        
        f_out = open(out, "w")
        
        for k, v in sorted(self.barcode_to_read.items()):
            q = sorted(v)
            f_out.write("%s\t%s\n" % (k, ",".join(q)))

        f_out.close()

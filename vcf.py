import sys
import pprint as pp

from helper import Helper


class VCF(Helper):

    def __init__(self, vcf, crossmap=None):
        self.vcf = vcf
        self.variants = {}
        self.variants_distribution = {}
        self.total = 0

        if crossmap:
            sys.path.append(crossmap)
            from CrossMap import read_chain_file
            from CrossMap import map_coordinates

            self.read_chain_file = read_chain_file
            self.map_coordinates = map_coordinates
    
    
    def format_af(self, i, AN):
        freq = float(i)/float(AN)
        if i == 0:
            return "%.2f" % (freq,)
        elif freq == 1.00:
            return "%.2f" % (freq,)
        elif freq >= 0.01:
            return "%.3f" % (freq,)
        else:
            return "%.3e" % (freq,)


    # TODO
    def cross_map(self, q_chr, q_start):
        
        q_end = q_start + 1

        # retrieve map file after loading chain file
        map_tree, target_chrom_sizes, src_chrom_sizes = self.read_chain_file(self.chain_file)

        # retrieve lifted over coordinates
        results = self.map_coordinates(mapTree, q_chr, q_start, q_end)

        if results:
            return (results[1][0], results[1][1])
        else:
            return results

    # TODO
    def cm_variants_from_ref_to_consensus_build(self, vcf, f_out):
        with open(vcf, "r") as f:
            for line in f:
                if line.startswith("#"):
                    f_out.write(line)
                else:
                    data = line.split("\t")
                    self.convert_to_consensus_build(data)

    # TODO
    def convert_to_consensus_build(self, data, ref, con):
        
        vcf = []
        chr = data[0]
        pos = int(data[1])
        var_id = data[2]
        ref_allele = data[3]
        alt_allele = data[4]
        qual = float(data[5])
        filter = data[6]
        info = data[7]
        format = data[8]
        genotypes = data[9:]

        # check whether reference has changed
        if ref_allele == ref[pos-1:pos-1+len(ref_allele)]:
            chr, pos = self.cross_map(data[0], int(data[1]))
            data[0] = chr
            data[1] = str(pos)
            return data  # change coordinates (if applicable), but no change in content
        else:
            pass

        # extract new ref base
        chr, pos = self.cross_map(data[0], int(data[1]))
        # convert alt bases

        # convert INFO

        # convert genotypes

    
    def get_variants(self, inferred_consensus, consensus_genomes):
        
        self.barcodes = sorted(consensus_genomes.keys())
        
        for barcode, cg in consensus_genomes.items():
            self.total += 1
            
            # search for all deletions
            del_position = None
            del_dict = {}
            for no, (x, y) in enumerate(zip(inferred_consensus,
                                            self.get_genome(cg))):
                if x != y:
                    if y == "D":
                        if del_position == None:
                            del_position = no + 1
                            del_indel = ["-"]
                        else:
                            del_indel.append("-")
                    else:
                        if del_position:
                            del_dict[del_position] = "".join(del_indel)
                        else:
                            del_position = None
            # last time
            if del_position:
                del_dict[del_position] = "".join(del_indel)

            counter = 0
            for no, (x, y) in enumerate(zip(inferred_consensus,
                                            self.get_genome(cg))):
                if x != y:
                    if x == "N":
                        continue
                    
                    if y == "D":
                        try:
                            pos = no + 1 - 1  # 0-based -> 1-based, and need to include ref/consensus base
                            var = inferred_consensus[no-1] + del_dict[no+1]
                        except KeyError:
                            continue
                    else:
                        pos = no + 1
                        var = y

                    # add variant to placeholder
                    try:
                        self.variants[pos]
                    except KeyError:
                        self.variants[pos] = {}

                    try:
                        self.variants[pos][var].append(self.barcodes.index(barcode))
                    except KeyError:
                        self.variants[pos][var] = [self.barcodes.index(barcode)]
                    
                    if y != "N":
                        counter += 1
                        
            try:
                self.variants_distribution[counter] += 1
            except:
                self.variants_distribution[counter] = 1
    
    
    def output_variants_distribution(self, f_out):
        
        f_out.write("##DistributionOfNumOfVariants\n")
        f_out.write("##" + "\t".join(["num_of_consensus_variants", "counts"]) + "\n")
        for num_of_var, counts in sorted(self.variants_distribution.items()):
            f_out.write("%s\t%s\n" % (num_of_var, counts))

    
    def output_vcf(self, inferred_consensus):
        #TODO variant quality
        #TODO genotype quality
        
        header = ["CHROM", "POS",
                  "ID", "REF",
                  "ALT", "QUAL",
                  "FILTER", "INFO",
                  "FORMAT"]
        header.extend(self.barcodes)
        f_out = open(self.vcf, "w")
        
        f_out.write("##fileformat=VCFv4.1\n")
        f_out.write("#%s\n" % ("\t".join(header)))

        for k, v in sorted(self.variants.items()):

            if ("N" in v.keys()) and (len(v.keys()) == 1):
                continue
            
            chr = "chrReference"
            pos = str(k)
            rsid = "."
            ref_base = inferred_consensus[k-1]
            alt_index = sorted([ alt for alt in v.keys() if alt != "N" ])
            alt_base = ",".join(alt_index)
            qual = "."
            filter = "."
            format = "GT"
            genotypes = [ 0 for i in range(len(self.barcodes)) ]
            info = [ "" for i in range(len(alt_index)) ]

            AC = [ 0 for i in alt_index ]
            AF = [ 0.0 for i in alt_index ]
            AN = self.total

            # work on correcting any deletion variants
            is_deletion = False
            no_of_deletions = []
            for var in v.keys():
                if "-" in var:
                    is_deletion = True
                    no_of_deletions.append(var.count("-"))
            
            if is_deletion:
                new_mapping = {}
                max_no_of_deletions = max(no_of_deletions)
                ref_base = inferred_consensus[k-1:k-1+max_no_of_deletions+1]
                temp = {}
                temp_index = []
                for var in v.keys():
                    if var != "N":
                        if "-" in var:  # deals with deletions
                            left_over = max_no_of_deletions - var.count("-")
                            if left_over == 0:
                                new_alt = ref_base[0]
                            else:
                                new_alt = ref_base[0] + ref_base[-left_over:]
                        else:  # deals with snps and insertions
                            new_alt = var + ref_base[1:]

                        new_mapping[var] = new_alt
                        temp_index.append(new_alt)

                # change alt_index
                alt_index = sorted(temp_index)

                # change alt_base
                alt_base = ",".join(alt_index)
                
                # change v
                for k2, v2 in v.items():
                    if k2 != "N":
                        temp[new_mapping[k2]] = v2
                    else:
                        temp[k2] = v2

                v = temp
            
            # work on genotypes
            for base, bcs in sorted(v.items()):

                if base == "N":
                    for bc in bcs:
                        genotypes[bc] = "."

                    AN -= len(bcs)
                else:
                    for bc in bcs:
                        genotypes[bc] = (alt_index.index(base) + 1)

            # work on INFO field
            for base, bcs in sorted(v.items()):

                if base != "N":
                    AC[alt_index.index(base)] = len(bcs)
                    AF[alt_index.index(base)] = self.format_af(len(bcs), AN)

            infov = ";".join(["AC=%s" % (",".join(map(str,AC))),
                              "AF=%s" % (",".join(map(str,AF))),
                              "AN=%s" % (AN,)])
            gt = "\t".join(map(str, genotypes))
            f_out.write("%s\n" % ("\t".join([chr, pos,
                                             rsid, ref_base,
                                             alt_base, qual,
                                             filter,
                                             infov,
                                             format,
                                             gt])))
        f_out.close()

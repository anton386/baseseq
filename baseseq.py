import sys
import time
import pysam
import pprint as pp
import getopt

from helper import Helper
from barcode import BarCode
from consensus import Consensus
from vcf import VCF
from chain import Chain

from lib.soma.ref import Reference

class BaseSeq(Helper):

    def __init__(self, bam, barcodes=None, out=None, ref=None,
                 rewritten_bam=None,
                 consensus_reference=None,
                 consensus_genomes=None,
                 haplotype_distribution=None,
                 vcf=None,
                 chain=None,
                 crossmap=None):
        
        self.bam = bam
        self.barcodes = barcodes
        self.out = out
        self.ref = Reference(ref)

        self.rewritten_bam = rewritten_bam
        self.rewritten_sorted_bam = rewritten_bam.replace(".bam", ".sorted.bam") if rewritten_bam else None
        self.consensus_reference = consensus_reference
        self.consensus_genomes = consensus_genomes
        self.haplotype_distribution = haplotype_distribution
        self.vcf = vcf
        self.chain = chain
        self.crossmap = crossmap
    
    def get_barcodes(self):
        # simple approach - align, take soft-clipped, and use the arbitrary 20 bases
        # intermediate approach - use the seed and extend approach
        out = open(self.out, "w")

        self.bc = BarCode(self.bam)
        self.bc.simple_approach()

        for k, v in sorted(self.bc.barcode_to_read.items()):
            q = sorted(v)
            out.write("%s\t%s\n" % (k, ",".join(q)))

        out.close()

    def error_correction_barcodes(self):

        # start analysis
        self.bc = BarCode(self.bam)
        sys.stdout.write("[%s] Starting Error Correction Analysis\n" % (self.get_time(),))

        # load barcodes
        self.bc.load_barcodes(self.barcodes)
        sys.stdout.write("[%s] Loaded BarCodes\n" % (self.get_time(),))

        # cluster barcodes
        self.bc.cluster_barcodes()
        sys.stdout.write("[%s] Clustered BarCodes\n" % (self.get_time(),))

    
    def filter_barcodes(self, barcode, export="fastq"):
        list_of_ids = []
        with open(self.barcodes, "r") as f:
            for line in f:
                data = line.strip("\r\n").split("\t")
                if barcode == data[0]:
                    list_of_ids = data[1].split(",")
                    break

        self.bc = BarCode(self.bam)
        self.bc.filter_and_export(list_of_ids, self.out, export=export)


    def sort_and_rewrite_bam(self):
        
        self.bc = BarCode(self.bam)
        sys.stderr.write("[%s] Starting Sort and Rewrite BAM\n" % (self.get_time(),))
        
        self.bc.load_barcodes(self.barcodes)
        sys.stderr.write("[%s] Loaded BarCodes\n" % (self.get_time(),))

        self.bc.sort_and_rewrite_bam(self.rewritten_bam)
        pysam.sort("-n", self.rewritten_bam, self.rewritten_sorted_bam.replace(".bam", ""))
        sys.stderr.write("[%s] Sort and Rewrite BAM\n" % (self.get_time(),))

    def split_bam_by_barcode(self):

        self.bc = BarCode(self.bam)
        sys.stderr.write("[%s] Starting procedure to split BAM by barcode\n" % (self.get_time(),))

        self.bc.split_bam_into_barcodes(self.out)
        sys.stderr.write("[%s] Finished splitting BAM by barcode id\n" % (self.get_time(),))
        

    def assemble_consensus_genomes(self):
        
        # build consensus
        self.consensus = Consensus(self.rewritten_sorted_bam, self.ref)
        sys.stderr.write("[%s] Starting Consensus Building\n" % (self.get_time(),))
        
        self.consensus.build()
        sys.stderr.write("[%s] Built and Calculated Consensus\n" % (self.get_time(),))
        
        self.consensus.infer_consensus(self.consensus_reference)
        sys.stderr.write("[%s] Inferred Consensus\n" % (self.get_time(),))

        self.consensus.output_consensus_genomes(self.consensus_genomes)
        sys.stderr.write("[%s] Output Consensus Genomes\n" % (self.get_time(),))

        self.consensus.output_haplotype_distribution(self.haplotype_distribution)
        sys.stderr.write("[%s] Output Haplotype Distribution\n" % (self.get_time(),))
        
        self.ovcf = VCF(self.vcf, crossmap=self.crossmap)
        #self.ovcf.get_variants(self.consensus.inferred_consensus,
        #                       self.consensus.consensus_genomes)
        #self.ovcf.output_vcf(self.consensus.inferred_consensus)
        self.ovcf.get_variants(self.ref.sequence,
                               self.consensus.consensus_genomes,
                               self.consensus.consensus_indels)
        self.ovcf.output_vcf(self.ref.sequence)
        sys.stderr.write("[%s] Output VCF\n" % (self.get_time(),))

        self.summary_statistics()
        sys.stderr.write("[%s] Output Summary Statistics\n" % (self.get_time(),))

        self.ochain = Chain(self.chain)
        self.ochain.output_chain(self.ref,
                                 self.consensus.inferred_consensus,
                                 self.consensus.inferred_fconsensus,
                                 self.consensus.inferred_consensus_indels)
        sys.stderr.write("[%s] Output Chain File\n" % (self.get_time(),))
    
    
    def summary_statistics(self):
        # coverage per genome
        # variants per genome
        # estimate PCR and sequencing errors
        # barcode distribution
        f_out = open(self.out, "w")

        self.bc = BarCode(self.bam)  #TEMP
        self.bc.load_barcodes(self.barcodes)  #TEMP
        
        self.consensus.output_consensus_coverage(f_out)
        self.ovcf.output_variants_distribution(f_out)
        self.bc.output_reads_in_barcode_distribution(f_out)
        
        f_out.close()

    def run(self):

        # Phase 1 - Detection of BarCode
        self.bc = BarCode(self.bam)
        sys.stderr.write("[%s] Starting BarCode Analysis \n" % (self.get_time(),))
        
        self.bc.simple_approach()
        sys.stderr.write("[%s] Analyzed BarCodes \n" % (self.get_time(),))
        
        self.bc.write_barcodes(self.barcodes)
        sys.stderr.write("[%s] Wrote BarCodes\n" % (self.get_time(),))

        # Phase 2 - Rewrite BAM
        sys.stderr.write("[%s] Starting Sort and Rewrite BAM\n" % (self.get_time(),))
        
        self.bc.load_barcodes(self.barcodes)
        sys.stderr.write("[%s] Loaded BarCodes\n" % (self.get_time(),))

        self.bc.bam.reset()
        self.bc.sort_and_rewrite_bam(self.rewritten_bam)
        pysam.sort("-n", self.rewritten_bam, self.rewritten_sorted_bam.replace(".bam", ""))
        sys.stderr.write("[%s] Sort and Rewrite BAM\n" % (self.get_time(),))
        
        # Phase 3 - Build Consensus
        self.consensus = Consensus(self.rewritten_sorted_bam, self.ref)
        sys.stderr.write("[%s] Starting Consensus Building\n" % (self.get_time(),))
        
        self.consensus.build()
        sys.stderr.write("[%s] Built and Calculated Consensus\n" % (self.get_time(),))
        
        self.consensus.infer_consensus(self.consensus_reference)
        sys.stderr.write("[%s] Inferred Consensus\n" % (self.get_time(),))

        # Phase 4 - Call Variants and Haplotypes
        self.consensus.output_consensus_genomes(self.consensus_genomes)
        sys.stderr.write("[%s] Output Consensus Genomes\n" % (self.get_time(),))

        self.consensus.output_haplotype_distribution(self.haplotype_distribution)
        sys.stderr.write("[%s] Output Haplotype Distribution\n" % (self.get_time(),))
        
        self.ovcf = VCF(self.vcf, crossmap=self.crossmap)
        self.ovcf.get_variants(self.ref.sequence,
                               self.consensus.consensus_genomes,
                               self.consensus.consensus_indels)
        self.ovcf.output_vcf(self.ref.sequence)
        sys.stderr.write("[%s] Output VCF\n" % (self.get_time(),))

        # Phase 5 - Summary Statistics and Chain Files
        f_out = open(self.out, "w")
        self.consensus.output_consensus_coverage(f_out)
        self.ovcf.output_variants_distribution(f_out)
        self.bc.output_reads_in_barcode_distribution(f_out)
        f_out.close()
        sys.stderr.write("[%s] Output Summary Statistics\n" % (self.get_time(),))

        self.ochain = Chain(self.chain)
        self.ochain.output_chain(self.ref,
                                 self.consensus.inferred_consensus,
                                 self.consensus.inferred_fconsensus,
                                 self.consensus.inferred_consensus_indels)
        sys.stderr.write("[%s] Output Chain File\n" % (self.get_time(),))
        
    
    def assemble_genomes(self):
        pass

    
    def assemble_genomes_from_fastq(self):
        pass


def usage():
    print """
python baseseq.py baseseq <options>
  -b --bam <non-sorted bam>
  -r --ref <reference>
  -a --barcodes <barcodes>
  -p --out-prefix <prefix of output or out directory>
  -c --crossmap-bin </path/to/dir>
  -h --help
"""

if __name__ == "__main__":

    method = sys.argv[1]

    if method == "get_barcodes":
        bam = sys.argv[2]
        out = sys.argv[3]
        
        bs = BaseSeq(bam, out=out)
        bs.get_barcodes()
    
    elif method == "filter_barcodes":
        bam = sys.argv[2]
        barcodes = sys.argv[3]
        barcode = sys.argv[6]
        out = sys.argv[4]
        export = sys.argv[5]
        
        bs = BaseSeq(bam, barcodes=barcodes, out=out)
        bs.filter_barcodes(barcode, export=export)

    elif method == "error_correction_barcodes":
        bam = sys.argv[2]
        barcodes = sys.argv[3]
        
        bs = BaseSeq(bam, barcodes=barcodes)
        bs.error_correction_barcodes()

    elif method == "consensus":
        bam = sys.argv[2]
        rewritten_bam = sys.argv[3]
        barcodes = sys.argv[4]
        ref = sys.argv[5]
        consensus_reference = sys.argv[6]
        consensus_genomes = sys.argv[7]
        haplotype_distribution = sys.argv[8]
        vcf = sys.argv[9]
        out = sys.argv[10]
        chain = sys.argv[11]

        bs = BaseSeq(bam, rewritten_bam=rewritten_bam,
                     barcodes=barcodes, ref=ref,
                     consensus_reference=consensus_reference,
                     consensus_genomes=consensus_genomes,
                     haplotype_distribution=haplotype_distribution,
                     vcf=vcf, out=out, chain=chain)
        bs.assemble_consensus_genomes()

    elif method == "sort_and_rewrite_bam":
        bam = sys.argv[2]
        barcodes = sys.argv[3]
        rewritten_bam = sys.argv[4]

        bs = BaseSeq(bam, barcodes=barcodes, rewritten_bam=rewritten_bam)
        bs.sort_and_rewrite_bam()

    elif method == "split_bam_by_barcode":
        bam = sys.argv[2]
        out = sys.argv[3]
        
        bs = BaseSeq(bam, out=out)
        bs.split_bam_by_barcode()

    elif method == "baseseq":

        try:
            short_options = "b:r:a:p:c:h"
            long_options = ["bam", "ref", "barcodes",
                            "out-prefix",
                            "crossmap-bin", "help"]
            opts, args = getopt.getopt(sys.argv[2:], short_options, long_options)
        except getopt.GetoptError as err:
            print str(err)
            usage()
            sys.exit(99)
        
        bam = None
        barcodes = None
        ref = None

        rewritten_bam = None
        consensus_reference = None
        consensus_genomes = None
        haplotype_distribution = None
        vcf = None
        out = None
        chain = None
        crossmap = None

        for o, a in opts:
            if o in ("-b", "--bam"):
                bam = a
            elif o in ("-r", "--ref"):
                ref = a
            elif o in ("-a", "--barcodes"):
                barcodes = a
            elif o in ("-c", "--crossmap-bin"):
                crossmap = a
            elif o in ("-p", "--out-prefix"):
                rewritten_bam = a + ".rewritten.bam"
                consensus_reference = a + ".consensus.fa"
                consensus_genomes = a + ".cg.fa"
                haplotype_distribution = a + ".freq"
                vcf = a + ".vcf"
                out = a + ".out"
                chain = a + ".chain"
            elif o in ("-h", "--help"):
                usage()
                sys.exit(99)

        bs = BaseSeq(bam, rewritten_bam=rewritten_bam,
                     barcodes=barcodes, ref=ref,
                     consensus_reference=consensus_reference,
                     consensus_genomes=consensus_genomes,
                     haplotype_distribution=haplotype_distribution,
                     vcf=vcf, out=out, chain=chain,
                     crossmap=crossmap)
        bs.run()

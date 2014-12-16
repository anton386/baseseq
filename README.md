# baseseq - Tool to generate haplotypes from BAsE-Seq libraries

## Acknowledgement
Hong et. al. - BAsE-Seq: a method for obtaining long viral haplotypes from short sequence reads

## Summary
This is an unofficial repository of scripts used to analyze viral HBV reads generated from the BAsE-Seq protocol.

## Prerequisite Dependencies

Python Libraries
  * pysam
  * CrossMap

## Usage

1. Run Analysis
```
python2.7 baseseq.py baseseq <options>
  -b --bam <bam>
  -r --ref <reference>
  -a --barcodes <barcodes>
  -p --out-prefix </path/to/prefix>
  -c --crossmap-bin </path/to/dir>
  -h --help
```

## Output Files
  * vcf (.vcf)
  * consensus genome (.consensus.fa)
  * demultiplexed consensus for individual genomes (.cg.fa)
  * haplotype frequency (.freq)
  * summary statistics (.out)
  * UCSC chain file between reference and consensus genome (.chain)
  * barcodes (.txt)

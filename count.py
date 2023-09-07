from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from Bio import SeqIO
from rich.console import Console
from rich.table import Table
from argparse import ArgumentParser
from itertools import islice
import sys
import gzip
import zlib


def read_fasta(fasta_file: str) -> dict:
    barcodes = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        length = len(record.seq)
        if length not in barcodes:
            barcodes[length] = set()
        barcodes[length].add(str(record.seq))
    return barcodes

def read_paired_fastq(fastq1_file: str, fastq2_file: str, chunk_size: int = 40000):
    paired_reads1 = []
    paired_reads2 = []
    with gzip.open(fastq1_file, 'rt') if fastq1_file.endswith('.gz') else open(fastq1_file, 'r') as f1, \
         gzip.open(fastq2_file, 'rt') if fastq2_file.endswith('.gz') else open(fastq2_file, 'r') as f2:
             
        while True:
            chunk1 = list(islice(f1, chunk_size))
            chunk2 = list(islice(f2, chunk_size))
            if not chunk1 or not chunk2:
                break
            
            paired_reads1.extend(chunk1[i+1].strip() for i in range(0, len(chunk1), 4))
            paired_reads2.extend(chunk2[i+1].strip() for i in range(0, len(chunk2), 4))
    
    return paired_reads1, paired_reads2

def process_chunk(args: tuple) -> Counter:
    chunk, barcodes_by_length = args
    counts = Counter()
    for rec1, rec2 in chunk:
        rec2_rev_comp = str(Seq(rec2).reverse_complement())
        for barcode_length in barcodes_by_length.keys():
            read_kmers_1 = {rec1[i:i + barcode_length] for i in range(len(rec1) - barcode_length + 1)}
            read_kmers_2 = {rec2_rev_comp[i:i + barcode_length] for i in range(len(rec2_rev_comp) - barcode_length + 1)}
            intersection = read_kmers_1 & read_kmers_2 & barcodes_by_length[barcode_length]
            counts.update(intersection)
    return counts

if __name__ == "__main__":
    parser = ArgumentParser(description='Process some barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input fasta file.')
    parser.add_argument('fastq1', type=str, help='First fastq file.')
    parser.add_argument('fastq2', type=str, help='Second fastq file.')
    args = parser.parse_args()

    # Initialize Rich console
    console = Console(stderr=True, highlight=False)
    console.rule("[bold red]Starting the program...")

    try:
        barcodes = read_fasta(args.fasta_file)
    except FileNotFoundError:
        console.print("File not found.")
        sys.exit(1)

    # Create a table for k-mer summary
    table = Table(show_header=True, header_style="bold green")
    table.add_column("K-mer Length", style="dim", width=12)
    table.add_column("Number of K-mers", style="dim", width=18)

    for length, kmer_set in barcodes.items():
        table.add_row(str(length), str(len(kmer_set)))

    console.rule("[bold red]K-mer Summary")
    console.print(table)

    # Read paired fastq files in chunks
    console.rule("[bold red]Reading fastq files in chunks...")
    chunk_size = 40000  # 10,000 reads x 4 lines per read
    reads1, reads2 = read_paired_fastq(args.fastq1, args.fastq2, chunk_size)

    all_reads = list(zip(reads1, reads2))

    num_cores = cpu_count()
    pool = Pool(processes=num_cores)

    chunk_size = len(all_reads) // num_cores
    reads_list = [(all_reads[i:i + chunk_size], barcodes) for i in range(0, len(all_reads), chunk_size)]

    results = pool.map(process_chunk, reads_list)

    final_counts = Counter()
    for res in results:
        final_counts.update(res)

    console.rule("[bold red]Final Barcode Counts")
    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")

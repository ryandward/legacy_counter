from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table
from argparse import ArgumentParser
import sys
import gzip

def read_fasta(fasta_file: str) -> dict:
    barcodes = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        length = len(record.seq)
        if length not in barcodes:
            barcodes[length] = set()
        barcodes[length].add(str(record.seq))
    return barcodes

def read_fastq(fastq_file: str) -> list:
    with gzip.open(fastq_file, 'rt') if fastq_file.endswith('.gz') else open(fastq_file, 'r') as f:
        return [str(record.seq) for record in SeqIO.parse(f, "fastq")]

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

    try:
        barcodes = read_fasta(args.fasta_file)
        reads1 = read_fastq(args.fastq1)
        reads2 = read_fastq(args.fastq2)
    except FileNotFoundError:
        print("File not found.")
        sys.exit(1)

    # Initialize Rich console
    console = Console(stderr=True, highlight=False)

    # Create a table
    table = Table(show_header=True, header_style="bold green")
    table.add_column("K-mer Length", style="dim", width=12)
    table.add_column("Number of K-mers", style="dim", width=18)

    for length, kmer_set in barcodes.items():
        table.add_row(str(length), str(len(kmer_set)))

    console.rule("[bold red]K-mer Summary")
    console.print(table)
    console.rule("[bold red]")

    all_reads = list(zip(reads1, reads2))

    num_cores = cpu_count()
    pool = Pool(processes=num_cores)

    chunk_size = len(all_reads) // num_cores
    reads_list = [(all_reads[i:i + chunk_size], barcodes) for i in range(0, len(all_reads), chunk_size)]

    results = pool.map(process_chunk, reads_list)

    final_counts = Counter()
    for res in results:
        final_counts.update(res)

    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")

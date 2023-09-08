from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from rich.console import Console
from argparse import ArgumentParser
import time
import gzip

def read_fasta(fasta_file: str) -> dict:
    with console.status("[bold green]Reading FASTA File..."):
        barcodes = {}
        open_func = gzip.open if fasta_file.endswith('.gz') else open
        with open_func(fasta_file, 'rt') as f:
            for line in f:
                if line.startswith(">"):
                    continue
                seq = line.strip()
                length = len(seq)
                if length not in barcodes:
                    barcodes[length] = set()
                barcodes[length].add(seq)
    return barcodes

def read_paired_fastq(fastq1_file: str, fastq2_file: str, num_threads: int):
    with console.status("[bold green]Reading FASTQ Files..."):
        paired_reads1 = []
        paired_reads2 = []
        
        open_func = gzip.open if fastq1_file.endswith('.gz') else open
        with open_func(fastq1_file, 'rt') as f1, open_func(fastq2_file, 'rt') as f2:
            while True:
                next(f1, None)  # Skip header
                next(f2, None)
                seq1_line = next(f1, None)
                seq2_line = next(f2, None)
                next(f1, None)  # Skip strand
                next(f2, None)
                next(f1, None)  # Skip quality
                next(f2, None)
                
                if seq1_line is None or seq2_line is None:
                    break
                
                paired_reads1.append(seq1_line.strip())
                paired_reads2.append(seq2_line.strip())

        chunk_size = len(paired_reads1) // num_threads
    return [(paired_reads1[i:i+chunk_size], paired_reads2[i:i+chunk_size], barcodes) for i in range(0, len(paired_reads1), chunk_size)]

def process_chunk(chunk: tuple) -> Counter:
    with console.status("[bold green]Processing Chunks..."):
        reads1, reads2, barcodes = chunk
        counts = Counter()
        for rec1, rec2 in zip(reads1, reads2):
            rec2_rev_comp = str(Seq(rec2).reverse_complement())
            for barcode_length in barcodes.keys():
                read_kmers_1 = {rec1[i:i + barcode_length] for i in range(len(rec1) - barcode_length + 1)}
                read_kmers_2 = {rec2_rev_comp[i:i + barcode_length] for i in range(len(rec2_rev_comp) - barcode_length + 1)}
                intersection = read_kmers_1 & read_kmers_2 & barcodes[barcode_length]
                counts.update(intersection)
    return counts

if __name__ == "__main__":
    parser = ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file. Can be gzipped.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file. Can be gzipped.')
    parser.add_argument('fastq2', type=str, help='Second FASTQ file. Can be gzipped.')
    args = parser.parse_args()

    console = Console(stderr=True, highlight=False)
    console.rule("[bold red]Starting the Program")

    start_time = time.time()
    barcodes = read_fasta(args.fasta_file)
    end_time = time.time()
    console.print(f"Time Taken to Read FASTA File: {end_time - start_time} Seconds")

    num_threads = cpu_count()

    start_time = time.time()
    chunks = read_paired_fastq(args.fastq1, args.fastq2, num_threads)
    end_time = time.time()
    console.print(f"Time Taken to Read FASTQ Files: {end_time - start_time} Seconds")

    start_time = time.time()
    pool = Pool(num_threads)
    results = pool.map(process_chunk, chunks)
    end_time = time.time()
    console.print(f"Time Taken to Process Chunks: {end_time - start_time} Seconds")

    final_counts = Counter()
    for res in results:
        final_counts.update(res)

    console.rule("[bold red]Final Barcode Counts")
    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")

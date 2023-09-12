from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from rich.console import Console
import argparse
import time
import gzip

def read_fasta(fasta_file):
    barcodes = set()
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            if not line.startswith(">"):
                barcodes.add(line.strip())
    return barcodes

import subprocess

def fastq_reader(proc_stdout):
    while True:
        next(proc_stdout, None)
        seq_line = next(proc_stdout, None)
        next(proc_stdout, None)
        next(proc_stdout, None)
        if seq_line is None:
            break
        yield seq_line.decode().strip()

def read_paired_fastq(fastq1_file, fastq2_file, num_threads):
    paired_reads1, paired_reads2 = [], []
    proc1 = proc2 = None
    if fastq1_file.endswith('.gz'):
        proc1 = subprocess.Popen(["pigz", "-dc", fastq1_file], stdout=subprocess.PIPE)
        proc2 = subprocess.Popen(["pigz", "-dc", fastq2_file], stdout=subprocess.PIPE)
        f1, f2 = proc1.stdout, proc2.stdout
    else:
        f1 = open(fastq1_file, 'rt')
        f2 = open(fastq2_file, 'rt')

    for seq1, seq2 in zip(fastq_reader(f1), fastq_reader(f2)):
        paired_reads1.append(seq1)
        paired_reads2.append(seq2)

    if proc1 and proc2:
        proc1.wait()
        proc2.wait()

    chunk_size = len(paired_reads1) // num_threads
    return [(paired_reads1[i:i+chunk_size], paired_reads2[i:i+chunk_size]) for i in range(0, len(paired_reads1), chunk_size)]



def process_chunk(chunk, barcodes, barcode_start1, barcode_start2, barcode_length):
    reads1, reads2 = chunk
    counts = Counter()
    for rec1, rec2 in zip(reads1, reads2):
        barcode1 = rec1[barcode_start1:barcode_start1 + barcode_length]
        barcode2 = str(Seq(rec2[barcode_start2:barcode_start2 + barcode_length]).reverse_complement())
        if barcode1 in barcodes and barcode2 in barcodes:
            counts[barcode1] += 1
    return counts

def find_start_positions(reads, barcodes, barcode_length, is_read2=False):
    for read in reads:
        for i in range(len(read) - barcode_length + 1):
            kmer = read[i:i+barcode_length]
            if is_read2:
                kmer = str(Seq(kmer).reverse_complement())
            if kmer in barcodes:
                return i

def main(args):
    num_threads = cpu_count()
    console = Console(stderr=True, highlight=False)
    console.rule("[bold red]Starting the Program")

    with console.status("[bold green]Reading FASTA File..."):
        start_time = time.time()
        barcodes = read_fasta(args.fasta_file)
        console.print(f"Time Taken to Read FASTA File: {time.time() - start_time} Seconds")

    with console.status("[bold green]Reading FASTQ Files and Finding Barcode Starts..."):
        start_time = time.time()
        chunks = read_paired_fastq(args.fastq1, args.fastq2, num_threads)
        sample1, sample2 = chunks[0]
        barcode_length = len(next(iter(barcodes)))
        barcode_start1 = find_start_positions(sample1, barcodes, barcode_length)
        barcode_start2 = find_start_positions(sample2, barcodes, barcode_length, is_read2=True)
        console.print(f"Time Taken: {time.time() - start_time} Seconds")

    with console.status("[bold green]Processing Chunks..."):
        start_time = time.time()
        pool = Pool(num_threads)
        results = pool.starmap(process_chunk, [(chunk, barcodes, barcode_start1, barcode_start2, barcode_length) for chunk in chunks])
        console.print(f"Time Taken to Process Chunks: {time.time() - start_time} Seconds")

    final_counts = Counter()
    for res in results:
        final_counts.update(res)

    total_reads = sum(len(chunk[0]) for chunk in chunks)
    num_barcodes_seen = len(final_counts)
    num_reads_with_barcode = sum(final_counts.values())
    
    console.rule("[bold red]Summary")
    console.print(f"Number of Barcodes: [bold]{len(barcodes)}[/bold]")
    console.print(f"Number of Reads: [bold]{total_reads}[/bold]")
    console.print(f"Number of Barcodes Seen: [bold]{num_barcodes_seen}[/bold]")
    console.print(f"Number of Reads with a Barcode: [bold]{num_reads_with_barcode}[/bold]")
    
    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file.')
    parser.add_argument('fastq2', type=str, help='Second FASTQ file.')
    args = parser.parse_args()

    main(args)

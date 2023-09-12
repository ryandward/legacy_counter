from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from rich.console import Console
import argparse
import time
import gzip
from datetime import timedelta

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

def determine_forward_read(sample1, sample2, barcodes):
    count1, count2 = 0, 0
    for seq1, seq2 in zip(sample1[:100], sample2[:100]):
        if any(barcode in seq1 for barcode in barcodes):
            count1 += 1
        if any(barcode in seq2 for barcode in barcodes):
            count2 += 1
    return count2 > count1

def process_chunk(chunk, barcodes, barcode_start1, barcode_start2, barcode_length):
    reads1, reads2 = chunk
    counts = Counter()
    unexpected_sequences = Counter()
    for rec1, rec2 in zip(reads1, reads2):
        barcode1 = rec1[barcode_start1:barcode_start1 + barcode_length]
        barcode2 = str(Seq(rec2[barcode_start2:barcode_start2 + barcode_length]).reverse_complement())
        if barcode1 in barcodes and barcode2 in barcodes:
            counts[barcode1] += 1
        else:
            unexpected_sequences[barcode1] += 1
    return counts, unexpected_sequences

def find_start_positions(reads, barcodes, barcode_length, is_read2=False):
    for read in reads:
        for i in range(len(read) - barcode_length + 1):
            kmer = read[i:i+barcode_length]
            if is_read2:
                kmer = str(Seq(kmer).reverse_complement())
            if kmer in barcodes:
                return i

from datetime import datetime
from rich.table import Table
from collections import Counter
from multiprocessing import Pool, cpu_count
from rich.console import Console

def main(args):
    num_threads = cpu_count()
    console = Console(stderr=True, highlight=False)
    console.rule("[bold red]Initializing Barcode Counting Operation")
    console.print(f"Program started at [bold]{datetime.now()}[/bold]")

    timing_table = Table(show_header=True, header_style="bold magenta")
    timing_table.add_column("Step", style="dim", width=50)
    timing_table.add_column("Time", justify="right")

    stats_table = Table(show_header=True, header_style="bold blue")
    stats_table.add_column("Statistic", style="dim", width=50)
    stats_table.add_column("Value", justify="right")

    start_global = datetime.now()

    start_time = datetime.now()
    with console.status("[bold green]Reading FASTA File..."):
        barcodes = read_fasta(args.fasta_file)
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Reading FASTA File", f"[bold]{last_step_time}")

    start_time = datetime.now()
    with console.status(f"[bold green]Reading FASTQ Files... (Last Step Took: {last_step_time})"):
        chunks = read_paired_fastq(args.fastq1, args.fastq2, num_threads)
        sample1, sample2 = chunks[0]
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Reading FASTQ Files", f"[bold]{last_step_time}")

    # Determine if forward direction needs to be swapped
    start_time = datetime.now()
    with console.status(f"[bold green]Determining Forward Direction..."):
        sample1, sample2 = chunks[0]
        need_swap = determine_forward_read(sample1, sample2, barcodes)  # Make sure the function name matches
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Determining Forward Direction", f"[bold]{last_step_time}")

    if need_swap:
        chunks = [(reads2, reads1) for reads1, reads2 in chunks]

    start_time = datetime.now()
    with console.status(f"[bold green]Determining Barcode Coordinates for Read 1... (Last Step Took: {last_step_time})"):
        barcode_length = len(next(iter(barcodes)))
        barcode_start1 = find_start_positions(sample1, barcodes, barcode_length)
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Determining Barcode Coordinates for Read 1", f"[bold]{last_step_time}")

    start_time = datetime.now()
    with console.status(f"[bold green]Determining Barcode Coordinates for Read 2... (Last Step Took: {last_step_time})"):
        barcode_start2 = find_start_positions(sample2, barcodes, barcode_length, is_read2=True)
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Determining Barcode Coordinates for Read 2", f"[bold]{last_step_time}")

    start_time = datetime.now()
    with console.status(f"[bold green]Processing Chunks... (Last Step Took: {last_step_time})"):
        pool = Pool(num_threads)
        results = pool.starmap(process_chunk, [(chunk, barcodes, barcode_start1, barcode_start2, barcode_length) for chunk in chunks])
    last_step_time = datetime.now() - start_time
    timing_table.add_row("Processing Chunks", f"[bold]{last_step_time}")

    total_time_taken = datetime.now() - start_global
    timing_table.add_row("Total Time Taken", f"[bold]{total_time_taken}")

    console.print(timing_table)

    final_counts = Counter()
    final_unexpected_sequences = Counter()
    for counts, unexpected_sequences in results:
        final_counts.update(counts)
        final_unexpected_sequences.update(unexpected_sequences)

    total_reads = sum(len(chunk[0]) for chunk in chunks)
    num_barcodes_seen = len(final_counts)
    num_reads_with_barcode = sum(final_counts.values())

    # Additional statistics
    most_frequent_barcode = final_counts.most_common(1)
    least_frequent_barcode = final_counts.most_common()[:-2:-1]
    fraction_reads_with_barcodes = num_reads_with_barcode / total_reads if total_reads > 0 else 0
    most_frequent_unexpected_seq = final_unexpected_sequences.most_common(1)

    # Add stats to the table
    stats_table.add_row("Barcode Start Location for Read 1", f"[bold]{barcode_start1}[/bold]")
    stats_table.add_row("Barcode Start Location for Read 2", f"[bold]{barcode_start2}[/bold]")
    stats_table.add_row("Number of Barcodes in Reference", f"[bold]{len(barcodes)}[/bold]")
    stats_table.add_row("Number of Unique Barcodes Seen", f"[bold]{num_barcodes_seen}[/bold]")
    stats_table.add_row("Total Number of Reads", f"[bold]{total_reads}[/bold]")
    stats_table.add_row("Number of Reads Containing a Barcode", f"[bold]{num_reads_with_barcode}[/bold]")
    stats_table.add_row("Most Frequent Barcode", f"[bold]{most_frequent_barcode}[/bold]")
    stats_table.add_row("Least Frequent Barcode", f"[bold]{least_frequent_barcode}[/bold]")
    stats_table.add_row("Fraction of Reads with Barcodes", f"[bold]{fraction_reads_with_barcodes:.2f}[/bold]")
    stats_table.add_row("Most Frequent Unexpected Sequence at Offset", f"[bold]{most_frequent_unexpected_seq}[/bold]")


    console.rule("[bold red]Summary")
    console.print(stats_table)
    console.print(f"Program finished at [bold]{datetime.now()}[/bold]")

    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file.')
    parser.add_argument('fastq2', type=str, help='Second FASTQ file.')
    args = parser.parse_args()

    main(args)

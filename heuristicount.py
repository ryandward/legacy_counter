import argparse
import subprocess
import gzip
import os
import sys

from datetime import datetime
from rich.table import Table
from rich.console import Console
from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from rich.console import Console


def read_fasta(fasta_file):
    barcodes = set()
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            if not line.startswith(">"):
                barcodes.add(line.strip())
    return barcodes

import subprocess

def fastq_reader(proc_stdout, should_decode):
    while True:
        next(proc_stdout, None)
        seq_line = next(proc_stdout, None)
        next(proc_stdout, None)
        next(proc_stdout, None)
        if seq_line is None:
            break
        yield seq_line.decode().strip() if should_decode else seq_line.strip()

def read_paired_fastq(fastq1_file, fastq2_file, num_threads):
    paired_reads1, paired_reads2 = [], []
    proc1 = proc2 = None

    if fastq1_file.endswith('.gz'):
        try:
            proc1 = subprocess.Popen(["pigz", "-dc", fastq1_file], stdout=subprocess.PIPE)
            proc2 = subprocess.Popen(["pigz", "-dc", fastq2_file], stdout=subprocess.PIPE)
            f1, f2 = proc1.stdout, proc2.stdout
        except FileNotFoundError:
            f1 = gzip.open(fastq1_file, 'rt')
            f2 = gzip.open(fastq2_file, 'rt')
    else:
        f1 = open(fastq1_file, 'rt')
        f2 = open(fastq2_file, 'rt')

    should_decode = bool(proc1 and proc2)

    for seq1, seq2 in zip(fastq_reader(f1, should_decode), fastq_reader(f2, should_decode)):
        paired_reads1.append(seq1)
        paired_reads2.append(seq2)

    if proc1 and proc2:
        proc1.wait()
        proc2.wait()

    chunk_size = len(paired_reads1) // num_threads
    return [(paired_reads1[i:i+chunk_size], paired_reads2[i:i+chunk_size]) for i in range(0, len(paired_reads1), chunk_size)]

import gzip  # Add this import for gzip support


def determine_forward_read(sample1, sample2, barcodes):
    count1, count2 = 0, 0
    barcode_length = len(next(iter(barcodes)))

    for seq1, seq2 in zip(sample1[:250], sample2[:250]):
        kmers1 = {seq1[i:i + barcode_length] for i in range(len(seq1) - barcode_length + 1)}
        kmers2 = {seq2[i:i + barcode_length] for i in range(len(seq2) - barcode_length + 1)}

        if barcodes & kmers1:
            count1 += 1
        if barcodes & kmers2:
            count2 += 1

    return count2 > count1

from collections import Counter

def find_ends(reads, start, length):
    lefts, rights = [], []

    for read in reads[:250]:
        left = read[start - 4: start]
        right = read[start + length: start + length + 4]
        lefts.append(left)
        rights.append(right)

    left = Counter(lefts).most_common(1)[0][0]
    right = Counter(rights).most_common(1)[0][0]
    return left, right

from collections import defaultdict

def process_chunk(chunk, barcodes, barcode_start1, barcode_start2, barcode_length):
    reads1, reads2 = chunk
    counts = defaultdict(int)
    unexpected_sequences = defaultdict(int)

    for rec1, rec2 in zip(reads1, reads2):
        candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
        candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]

        if candidate1 in barcodes and candidate2 and candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN")) == candidate1:
            counts[candidate1[4:-4]] += 1
        else:
            if candidate1 not in barcodes and candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN")) == candidate1:
                unexpected_sequences[candidate1] += 1
    return counts, unexpected_sequences

def find_start_positions(reads, barcodes, barcode_length, is_read2=False):
    offset_counts = Counter()
    for read in reads[:250]:
        for i in range(len(read) - barcode_length + 1):
            kmer = read[i:i+barcode_length]
            if is_read2:
                kmer = str(Seq(kmer).reverse_complement())
            if kmer in barcodes:
                offset_counts[i] += 1
    return offset_counts.most_common(1)[0][0]

def timed_action(console, table, action, description, start_global=None):
    start_time = datetime.now()
    time_since_global = datetime.now() - start_global if start_global else None
    desc = f"{description} (Began After: {time_since_global})" if time_since_global else description
    with console.status(f"[bold green]{desc}..."):
        action_result = action()
    last_step_time = datetime.now() - start_time
    table.add_row(description, f"[bold]{last_step_time}")
    return action_result

def main(args):
    num_threads = cpu_count()
    console = Console(stderr=True, highlight=False)
    console.rule("[bold red]Initializing Barcode Counting Operation")
    console.print(f"Program started at [bold]{datetime.now()}[/bold]")

    timing_table = Table(show_header=True, header_style="bold magenta")
    timing_table.add_column("Step", style="dim", width=50)
    timing_table.add_column("Time", justify="right")

    start_global = datetime.now()

    barcodes = timed_action(console, timing_table, lambda: read_fasta(args.fasta_file), "Reading FASTA File", start_global)
    barcode_length = len(next(iter(barcodes)))

    chunks = timed_action(console, timing_table, lambda: read_paired_fastq(args.fastq1, args.fastq2, num_threads), "Reading FASTQ Files", start_global)
    sample1, sample2 = chunks[0]

    need_swap = timed_action(console, timing_table, lambda: determine_forward_read(sample1, sample2, barcodes), "Finding Forward Direction", start_global)

    if need_swap:
        chunks = [(reads2, reads1) for reads1, reads2 in chunks]
        sample1, sample2 = chunks[0]

    barcode_start1 = timed_action(console, timing_table, lambda: find_start_positions(sample1, barcodes, barcode_length), "Finding Forward Coordinates", start_global)
    barcode_start2 = timed_action(console, timing_table, lambda: find_start_positions(sample2, barcodes, barcode_length, is_read2=True), "Finding Reverse Coordinates", start_global)

    left1, right1 = timed_action(console, timing_table, lambda: find_ends(sample1, barcode_start1, barcode_length), "Finding Forward Junctions", start_global)
    left2, right2 = timed_action(console, timing_table, lambda: find_ends(sample2, barcode_start2, barcode_length), "Finding Reverse Junctions", start_global)

    # Adjust coordinates
    barcode_start1 -= len(left1)
    barcode_start2 -= len(left2)

    # Update barcodes
    barcodes = {left1 + barcode + right1 for barcode in barcodes}
    barcode_length = len(next(iter(barcodes)))

    with Pool(num_threads) as pool:
        results = timed_action(console, timing_table, lambda: pool.starmap(process_chunk, [(chunk, barcodes, barcode_start1, barcode_start2, barcode_length) for chunk in chunks]), "Processing Chunks", start_global)

    final_counts = Counter()
    final_unexpected_sequences = Counter()

    for counts, unexpected_sequences in results:
        final_counts.update(counts)
        final_unexpected_sequences.update(unexpected_sequences)

    # This checks the final_unexpected_sequences and checks if they each begin with left1 and end with right1, then adds them to a new Counter called new_barcodes
    new_barcodes = Counter()
    for seq, count in final_unexpected_sequences.items():
        if seq.startswith(left1) and seq.endswith(right1):
            new_barcodes[seq[4:-4]+"*"] = count

    final_counts = final_counts + new_barcodes

    total_reads = sum(len(chunk[0]) for chunk in chunks)
    num_barcodes_seen = len(final_counts)
    num_reads_with_barcode = sum(final_counts.values())

    # Additional statistics
    
    most_frequent_barcode = final_counts.most_common(1)
    fraction_reads_with_barcodes = num_reads_with_barcode / total_reads if total_reads > 0 else 0
    most_frequent_unexpected_seq = final_unexpected_sequences.most_common(1)
    fastq1_filename = os.path.basename(args.fastq1) if not need_swap else os.path.basename(args.fastq2)
    fastq2_filename = os.path.basename(args.fastq2) if not need_swap else os.path.basename(args.fastq1)

    # Tables

    # Table for input and configuration
    input_config_table = Table(show_header=True, header_style="bold blue")
    input_config_table.add_column("Input & Configuration", style="dim", width=50)
    input_config_table.add_column("Value", justify="right")
    input_config_table.add_row("Forward Read File", f"[bold]{fastq1_filename}[/bold]")
    input_config_table.add_row("Reverse Read File", f"[bold]{fastq2_filename}[/bold]")
    input_config_table.add_row("Barcodes File", f"[bold]{os.path.basename(args.fasta_file)}[/bold]")
    input_config_table.add_row("Barcode Start Location for Forward", f"[bold]{barcode_start1}[/bold]")
    input_config_table.add_row("Barcode Start Location for Reverse", f"[bold]{barcode_start2}[/bold]")

    # Table for numeric statistics
    numeric_stats_table = Table(show_header=True, header_style="bold green")
    numeric_stats_table.add_column("Numeric Statistic", style="dim", width=50)
    numeric_stats_table.add_column("Value", justify="right")
    numeric_stats_table.add_row("Number of Barcodes in Reference", f"[bold]{len(barcodes)}[/bold]")
    numeric_stats_table.add_row("Number of Unique Barcodes Seen", f"[bold]{num_barcodes_seen}[/bold]")
    numeric_stats_table.add_row("Total Number of Reads", f"[bold]{total_reads}[/bold]")
    numeric_stats_table.add_row("Number of Reads Containing a Barcode", f"[bold]{num_reads_with_barcode}[/bold]")
    numeric_stats_table.add_row("Fraction of Reads with Barcodes", f"[bold]{fraction_reads_with_barcodes:.2f}[/bold]")

    # Table for sequence information
    sequence_info_table = Table(show_header=True, header_style="bold yellow")
    sequence_info_table.add_column("Sequence Information", style="dim", width=50)
    sequence_info_table.add_column("Value", justify="right")
    sequence_info_table.add_row("Most Frequent Barcode", f"[bold]{most_frequent_barcode}[/bold]")
    sequence_info_table.add_row("Most Frequent Unexpected Sequence", f"[bold]{new_barcodes.most_common(1)}[/bold]")

    # Finish Timing
    total_time_taken = datetime.now() - start_global
    timing_table.add_row("Total Time Taken", f"[bold]{total_time_taken}")

    # Print Tables
    console.print(input_config_table)
    console.print(timing_table)
    console.print(numeric_stats_table)
    console.print(sequence_info_table)

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
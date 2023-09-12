from multiprocessing import Pool, cpu_count
from collections import Counter
from Bio.Seq import Seq
from rich.console import Console
from argparse import ArgumentParser
import time
import gzip

console = Console(stderr=True, highlight=False)

def read_fasta(fasta_file: str) -> set:
    """Read barcodes from a FASTA file."""
    barcodes = set()
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq = line.strip()
            barcodes.add(seq)
    return barcodes

def fastq_reader(file):
    """Generator function to read records from a FASTQ file."""
    while True:
        next(file, None)  # Skip header
        seq_line = next(file, None)  # Read sequence line
        next(file, None)  # Skip strand
        next(file, None)  # Skip quality
        if seq_line is None:
            break
        yield seq_line.strip()

def read_paired_fastq(fastq1_file: str, fastq2_file: str, num_threads: int):
    """Read paired-end sequences from two FASTQ files and chunk them for parallel processing."""
    paired_reads1, paired_reads2 = [], []
    
    open_func = gzip.open if fastq1_file.endswith('.gz') else open
    with open_func(fastq1_file, 'rt') as f1, open_func(fastq2_file, 'rt') as f2:
        reader1, reader2 = fastq_reader(f1), fastq_reader(f2)
        for seq1, seq2 in zip(reader1, reader2):
            paired_reads1.append(seq1)
            paired_reads2.append(seq2)
    
    chunk_size = len(paired_reads1) // num_threads
    return [
        (paired_reads1[i:i+chunk_size], paired_reads2[i:i+chunk_size], barcodes)
        for i in range(0, len(paired_reads1), chunk_size)
    ]

def determine_forward_direction(chunk: tuple, sample_size: int = 100) -> str:
    """Determine which direction is forward based on the intersection of barcodes and reads."""
    reads1, reads2, barcodes = chunk
    barcode_length = len(next(iter(barcodes)))
    
    count1, count2 = 0, 0
    for i, (rec1, rec2) in enumerate(zip(reads1, reads2)):
        if i >= sample_size:
            break

        kmer_set1 = {rec1[j:j + barcode_length] for j in range(len(rec1) - barcode_length + 1)}
        kmer_set2 = {rec2[j:j + barcode_length] for j in range(len(rec2) - barcode_length + 1)}
        
        count1 += len(kmer_set1 & barcodes)
        count2 += len(kmer_set2 & barcodes)

    return "forward" if count1 > count2 else "reverse"

def process_chunk(chunk: tuple) -> Counter:
    """Process a chunk of paired-end reads to count barcode occurrences."""
    reads1, reads2, barcodes = chunk
    counts = Counter()
    barcode_length = len(next(iter(barcodes)))
    
    cache = {}  # Cache for storing counts of previously seen read pairs

    for rec1, rec2 in zip(reads1, reads2):
        rec_pair = (rec1, rec2)
        
        if rec_pair in cache:
            counts.update(cache[rec_pair])
            continue

        # Reverse complement the second read
        rec2_rev_comp = str(Seq(rec2).reverse_complement())
        
        # Generate k-mers from the reads
        kmer_set1 = {rec1[i:i + barcode_length] for i in range(len(rec1) - barcode_length + 1)}
        kmer_set2 = {rec2_rev_comp[i:i + barcode_length] for i in range(len(rec2_rev_comp) - barcode_length + 1)}
        
        # Early termination if no intersection with barcodes
        if not (kmer_set1 & barcodes) or not (kmer_set2 & barcodes):
            continue
        
        intersect_kmers = kmer_set1 & kmer_set2 & barcodes
        cache[rec_pair] = intersect_kmers

        for kmer in intersect_kmers:
            counts[kmer] += 1

    return counts

if __name__ == "__main__":
    parser = ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file.')
    parser.add_argument('fastq2', type=str, help='Second FASTQ file.')
    args = parser.parse_args()

    num_threads = cpu_count()
    console.rule("[bold red]Starting the Program")

    with console.status("[bold green]Reading FASTA File..."):
        start_time = time.time()
        barcodes = read_fasta(args.fasta_file)
        elapsed_time = round(time.time() - start_time, 2)
        console.print(f"Time Taken to Read FASTA: {elapsed_time}s. Number of Barcodes: [bold]{len(barcodes)}[/bold]")

    with console.status("[bold green]Reading FASTQ Files..."):
        start_time = time.time()
        chunks = read_paired_fastq(args.fastq1, args.fastq2, num_threads)
        elapsed_time = round(time.time() - start_time, 2)
        total_reads = sum(len(chunk[0]) for chunk in chunks)
        console.print(f"Time Taken to Read FASTQ: {elapsed_time}s. Number of Reads: [bold]{total_reads}[/bold]")

    with console.status("[bold green]Determining Forward Direction..."):
        start_time = time.time()
        sample_chunk = chunks[0]
        direction = determine_forward_direction(sample_chunk)
        elapsed_time = round(time.time() - start_time, 2)
        console.print(f"Time Taken to Determine Direction: {elapsed_time}s. Forward Direction: [bold]{direction}[/bold]")

    with console.status("[bold green]Processing Chunks..."):
        start_time = time.time()
        pool = Pool(num_threads)
        if direction == "forward":
            results = pool.map(process_chunk, chunks)
        else:
            reversed_chunks = [(reads2, reads1, barcodes) for reads1, reads2, barcodes in chunks]
            results = pool.map(process_chunk, reversed_chunks)
        elapsed_time = round(time.time() - start_time, 2)
        console.print(f"Time Taken to Process Chunks: {elapsed_time}s")

    final_counts = Counter()
    for res in results:
        final_counts.update(res)

    num_barcodes_seen = len(final_counts)
    num_reads_with_barcode = sum(final_counts.values())
    
    console.rule("[bold red]Summary")
    console.print(f"Number of Barcodes: [bold]{len(barcodes)}[/bold]")
    console.print(f"Number of Reads: [bold]{total_reads}[/bold]")
    console.print(f"Number of Barcodes Seen: [bold]{num_barcodes_seen}[/bold]")
    console.print(f"Number of Reads with a Barcode: [bold]{num_reads_with_barcode}[/bold]")

    console.rule("[bold red]Final Barcode Counts")
    for barcode, count in final_counts.items():
        print(f"{barcode}\t{count}")
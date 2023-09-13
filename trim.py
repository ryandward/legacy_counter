import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

def open_file(file_path):
    return gzip.open(file_path, 'rt') if file_path.endswith('.gz') else open(file_path, 'r')

def load_barcodes(fasta_file):
    barcodes = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        barcodes.add(str(record.seq))
    return barcodes

def orientation_check(fastq_file, barcodes, sample_size=20):
    forward_count = 0
    reverse_count = 0
    for idx, record in enumerate(SeqIO.parse(fastq_file, "fastq")):
        if idx >= sample_size:
            break
        sequence = str(record.seq)
        for barcode in barcodes:
            if barcode in sequence:
                forward_count += 1
            elif str(Seq(barcode).reverse_complement()) in sequence:
                reverse_count += 1
    return "reverse" if reverse_count > forward_count else "forward"

def find_flanking_sequences(fastq_file, barcodes, is_reverse, sample_size=20):
    flanking_sequences = {"upstream": [], "downstream": []}
    for idx, record in enumerate(SeqIO.parse(fastq_file, "fastq")):
        if idx >= sample_size:
            break
        sequence = str(record.seq) if is_reverse == "forward" else str(Seq(str(record.seq)).reverse_complement())
        for barcode in barcodes:
            start_idx = sequence.find(barcode)
            if start_idx == -1:
                continue
            end_idx = start_idx + len(barcode)
            upstream = sequence[max(0, start_idx - 5):start_idx]
            downstream = sequence[end_idx:end_idx + 5]
            flanking_sequences["upstream"].append(upstream)
            flanking_sequences["downstream"].append(downstream)

    common_upstream = Counter(flanking_sequences["upstream"]).most_common(1)[0][0]
    common_downstream = Counter(flanking_sequences["downstream"]).most_common(1)[0][0]
    return common_upstream, common_downstream

def trim_and_output(fastq_file, barcodes, upstream, downstream, is_reverse):
    for record in SeqIO.parse(fastq_file, "fastq"):
        sequence = str(record.seq) if is_reverse == "forward" else str(record.seq.reverse_complement())
        found = False
        for barcode in barcodes:
            full_pattern = upstream + barcode + downstream
            start_idx = sequence.find(full_pattern)
            if start_idx == -1:
                continue
            end_idx = start_idx + len(full_pattern)
            trimmed_start = start_idx + len(upstream)
            trimmed_end = end_idx - len(downstream)
            trimmed_record = record[trimmed_start:trimmed_end]
            trimmed_record.id = record.id
            trimmed_record.description = ""
            SeqIO.write(trimmed_record, sys.stdout, "fastq")
            found = True
            break

        if not found:
            print(f"@{record.id}")
            print()
            print("+")
            print()

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    fastq_file_path = sys.argv[2]

    barcodes = load_barcodes(fasta_file)
    is_reverse = orientation_check(open_file(fastq_file_path), barcodes)
    upstream, downstream = find_flanking_sequences(open_file(fastq_file_path), barcodes, is_reverse)
    trim_and_output(open_file(fastq_file_path), barcodes, upstream, downstream, is_reverse)

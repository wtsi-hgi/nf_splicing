import re
import subprocess

def pigz_open(path: str):
    """
    Open a gzip file using pigz for faster decompression
    Parameters:
        -- path: file path to the gzip file
    Returns:
        -- io_wrapper: a TextIOWrapper for reading the decompressed file
    """
    return subprocess.Popen(["pigz", "-dc", path], stdout = subprocess.PIPE)

def fastq_iter(handle):
    """
    FASTQ parser yielding (header, seq, qual)
    Parameters:
        -- handle: file handle for the FASTQ file
    Yields:
        -- (header, seq, qual): a tuple containing the header, sequence, and quality
    """
    while True:
        header = handle.readline()
        if not header:
            break
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not (seq and plus and qual):
            raise ValueError("Truncated FASTQ record")

        yield (header.rstrip("\n"), seq.rstrip("\n"), qual.rstrip("\n"))

def read_first_fasta_seq(fasta_path):
    """
    Read the first record of the fasta file
    Parameters:
        -- fasta_path: the path of fasta file
    Returns:
        -- str: the first sequence of the fasta file
    """
    seq_lines = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if seq_lines:
                    break
                continue
            seq_lines.append(line)
    return "".join(seq_lines)

def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.
    Parameters:
        -- seq: the DNA sequence to reverse complement
    Returns:
        -- str: the reverse complement of the sequence
    """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def hamming_distance(str1: str, str2: str) -> int:
    """
    Calculate the Hamming distance between two strings
    Parameters:
        -- str1: first string
        -- str2: second string
    Returns:
        -- int: the Hamming distance, or the maximum length if they differ in length
    """
    if len(str1) != len(str2):
        return max(len(str1), len(str2))
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

def match_approximate(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq allowing max_mismatches.
    Returns start index of match or -1 if not found.
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    k = len(pattern)
    for i in range(len(seq) - k + 1):
        window = seq[i:i+k]
        if hamming_distance(window, pattern) <= max_mismatches:
            return i
    return -1

def extract_sequence(seq: str, up_seq: str, down_seq: str, max_mismatches: int) -> str:
    """
    Extract substring from seq between approximate matches of up_seq and down_seq.
    Parameters:
        -- seq: the sequence to search in
        -- up_seq: upstream sequence to match
        -- down_seq: downstream sequence to match
        -- max_mismatches: maximum number of mismatches allowed for both matches
    Returns:
        -- str: the extracted sequence or an error message if not found
    """
    start_idx = match_approximate(seq, up_seq, max_mismatches)
    if start_idx == -1:
        return "upstream not found"
    start_idx += len(up_seq)

    end_idx = match_approximate(seq[start_idx:], down_seq, max_mismatches)
    if end_idx == -1:
        return "downstream not found"
    end_idx += start_idx

    return seq[start_idx:end_idx]

def check_barcode(barcode_seq: str, barcode_temp: str, max_mismatches: int):
    """
    Check if barcode_seq matches the barcode_temp allowing max_mismatches
    Parameters:
        -- barcode_seq: the sequence to check
        -- barcode_temp: the template sequence to match against
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- str or None: corrected barcode sequence if matches, else None
    """
    if len(barcode_seq) != len(barcode_temp):
        return None

    mismatch_count = 0
    barcode_corrected = []

    for s_char, p_char in zip(barcode_seq, barcode_temp):
        if p_char == 'N':
            barcode_corrected.append(s_char)
        elif s_char != p_char:
            mismatch_count += 1
            barcode_corrected.append(p_char)
            if mismatch_count > max_mismatches:
                return None
        else:
            barcode_corrected.append(s_char)

    return "".join(barcode_corrected)

def calc_softclip_lens(cigar: str) -> tuple[int, int]:
    first_softclip = 0
    last_softclip = 0

    match_start = re.match(r'^(\d+)S', cigar)
    if match_start:
        first_softclip = int(match_start.group(1))

    match_end = re.search(r'(\d+)S$', cigar)
    if match_end:
        last_softclip = int(match_end.group(1))

    return first_softclip, last_softclip
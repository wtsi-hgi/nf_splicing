import re
import numpy as np
import subprocess
import Levenshtein

def pigz_open(path: str):
    """
    Open a gzip file using pigz for faster decompression
    Parameters:
        -- path: file path to the gzip file
    Returns:
        -- io_wrapper: a TextIOWrapper for reading the decompressed file
    """
    return subprocess.Popen(["pigz", "-dc", path], stdout = subprocess.PIPE)

def fastq_iter_se(handle):
    """
    FASTQ single-end parser yielding (header, seq, qual)
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

def fastq_iter_pe(handle1, handle2):
    """
    FASTQ paired-end parser yielding ((header1, seq1, qual1), (header2, seq2, qual2))
    Parameters:
        -- handle1: file handle for read 1 FASTQ
        -- handle2: file handle for read 2 FASTQ
    Yields:
        -- ((header1, seq1, qual1), (header2, seq2, qual2))
    """
    while True:
        header1 = handle1.readline()
        header2 = handle2.readline()
        if not header1 or not header2:
            break

        seq1 = handle1.readline()
        seq2 = handle2.readline()
        plus1 = handle1.readline()
        plus2 = handle2.readline()
        qual1 = handle1.readline()
        qual2 = handle2.readline()

        # check for truncated records
        if not (seq1 and plus1 and qual1 and seq2 and plus2 and qual2):
            raise ValueError("Truncated FASTQ record in one of the pairs")

        yield ((header1.rstrip("\n"), seq1.rstrip("\n"), qual1.rstrip("\n")), 
               (header2.rstrip("\n"), seq2.rstrip("\n"), qual2.rstrip("\n")))

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

# def hamming_distance(str1: str, str2: str) -> int:
#     """
#     Calculate the Hamming distance between two strings
#     Parameters:
#         -- str1: first string
#         -- str2: second string
#     Returns:
#         -- int: the Hamming distance, or the maximum length if they differ in length
#     """
#     str1_bytes = str1.encode("ascii")
#     str2_bytes = str2.encode("ascii")
#     return sum(b1 != b2 for b1, b2 in zip(str1_bytes, str2_bytes))

# def match_approximate(seq: str, pattern: str, max_mismatches: int, distance: str) -> int:
#     """
#     Find approximate match of pattern in seq allowing max_mismatches.
#     Returns start index of match or -1 if not found.
#     Parameters:
#         -- seq: the sequence to search in
#         -- pattern: the pattern to match
#         -- max_mismatches: maximum number of mismatches allowed
#     Returns:
#         -- int: start index of the match or -1 if not found
#     """
#     k = len(pattern)
#     n = len(seq)

#     if k > n:
#         return -1

#     # Select distance function
#     if distance == "hamming":
#         dist_func = hamming_distance
#     elif distance == "levenshtein":
#         dist_func = Levenshtein.distance
#     else:
#         raise ValueError(f"Unknown distance metric: {distance}")

#     # Slide window
#     for i in range(n - k + 1):
#         window = seq[i:i+k]
#         if dist_func(window, pattern) <= max_mismatches:
#             return i

#     return -1

def match_hamming_numpy(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq by hamming distance allowing max_mismatches
    Parameters:
        -- seq: the target sequence
        -- pattern: the pattern to match
    Returns:
        -- int: the Hamming distance, or the maximum length if they differ in length
    """
    k = len(pattern)
    n = len(seq)
    if k > n:
        return -1

    # convert the sequence to np.uint8 arrays of ASCII codes
    seq_arr = np.frombuffer(seq.encode("ascii"), dtype = np.uint8)
    pat_arr = np.frombuffer(pattern.encode("ascii"), dtype = np.uint8)

    # sliding window: create a 2D view of seq_arr of shape (n-k+1, k)
    windows = np.lib.stride_tricks.sliding_window_view(seq_arr, window_shape = k)

    # calculate hamming distances vectorized
    mismatches = np.sum(windows != pat_arr, axis = 1)
    matches = np.where(mismatches <= max_mismatches)[0]
    return int(matches[0]) if matches.size > 0 else -1

def match_levenshtein(seq: str, pattern: str, max_mismatches: int) -> int:
    """
    Find approximate match of pattern in seq by levenshtein distance allowing max_mismatches
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    k = len(pattern)
    n = len(seq)
    if k > n:
        return -1

    for i in range(n - k + 1):
        window = seq[i:i+k]
        if Levenshtein.distance(window, pattern) <= max_mismatches:
            return i
    return -1

def match_approximate(seq: str, pattern: str, max_mismatches: int, distance: str) -> int:
    """
    Hybrid approximate match supporting Hamming and Levenshtein
    Parameters:
        -- seq: the sequence to search in
        -- pattern: the pattern to match
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- int: start index of the match or -1 if not found
    """
    if distance == "hamming":
        return match_hamming_numpy(seq, pattern, max_mismatches)
    elif distance == "levenshtein":
        return match_levenshtein(seq, pattern, max_mismatches)
    else:
        raise ValueError(f"Unknown distance metric: {distance}")

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
    start_idx = match_approximate(seq, up_seq, max_mismatches, "hamming")
    if start_idx == -1:
        return "upstream not found"
    start_idx += len(up_seq)

    end_idx = match_approximate(seq[start_idx:], down_seq, max_mismatches, "hamming")
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
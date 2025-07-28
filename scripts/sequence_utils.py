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

def check_barcode(barcode_seq: str, barcode_temp: str, max_mismatches: int) -> bool:
    """
    Check if barcode_seq matches the barcode_temp allowing max_mismatches
    Parameters:
        -- barcode_seq: the sequence to check
        -- barcode_temp: the template sequence to match against
        -- max_mismatches: maximum number of mismatches allowed
    Returns:
        -- bool: True if matches within allowed mismatches, False otherwise
    """
    if len(barcode_seq) != len(barcode_temp):
        return False

    mismatch_count = 0
    for s_char, p_char in zip(barcode_seq, barcode_temp):
        if p_char != 'N':
            if s_char != p_char:
                mismatch_count += 1
                if mismatch_count > max_mismatches:
                    return False

    return True
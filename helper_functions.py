import math

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    compared_positions = 0
    mismatches = 0

    # go through both sequences together, one position at a time
    for i in range(0, len(reference_sequence)):
        base_ref = reference_sequence[i]
        base_seq = distant_sequence[i]

        if base_ref == '-' or base_seq == '-':
            continue

        compared_positions = compared_positions + 1
        if base_ref != base_seq:
            mismatches = mismatches + 1

    # raw proportion of differences
    proportion_diff = mismatches / compared_positions
    # Jukes-Cantor correction
    distance_per_nuc = -3/4 * math.log(1 - 4/3 * proportion_diff)
    # convert to per genome
    distance_per_genome = distance_per_nuc * compared_positions
    return distance_per_genome


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    raise NotImplementedError()

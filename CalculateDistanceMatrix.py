import math
from typing import Dict, Tuple, List
import numpy as np


def calculate_p_distance(sq1: str, sq2: str) -> float:
    """
   Calculate the proportion of differing sites between two sequences.
   The program uses the direct analytical solution (the third option mentioned in the requirements
   - "להיות חכמים יותר"). Instead of using line search or grid search, it uses the Jukes-Cantor formula which gives us the
   maximum likelihood estimate directly:
    d = -(3/4)ln(1 - 4p/3)
    Where:
    p is the p-distance (proportion of different sites)
    d is the evolutionary distance in time units

    This is optimal because:

    Under the Jukes-Cantor model with rate matrix R as given (all substitutions have rate 1/3), this formula
     gives the maximum likelihood estimate of the evolutionary time
    It's computationally efficient - no need for iterative optimization
    It's exact - no approximation like in grid search

    The program handles the case where p is too large (> 0.75) by returning infinity, which is appropriate
     since such large distances suggest saturation of substitutions.
   """
    if len(sq1) != len(sq2):
        raise ValueError("Sequences must be of equal length")
    differences = sum(1 for a, b in zip(sq1, sq2) if a != b)

    return differences / len(sq2)


def jukes_cantor_distance(p: float) -> float:
    """
    Calculate Jukes-Cantor distance from p-distance.
    Using the formula: d = -3/4 * ln(1 - 4/3 * p)
    """
    # Handle edge cases
    if p >= 0.75:  # theoretical maximum for JC69
        return float('inf')

    try:
        distance = -0.75 * math.log(1 - (4.0 / 3.0) * p, math.e)
        return distance
    except ValueError:  # Handle numerical issues
        return float('inf')


def calculate_dist_matrix(sequences: Dict[str, str]) -> Tuple[np.ndarray, List[str]]:
    """
    Calculate the distance matrix using Jukes-Cantor model.
    Returns the matrix and the list of sequence IDs in order.
    sequences = {id: sequence}
    """
    sequences_id = list(sequences.keys())
    num_of_seq = len(sequences_id)

    #initialize dist matrix
    dist_matrix = np.zeros((num_of_seq, num_of_seq))

    #find dist for each 2 sequences
    for i in range(num_of_seq):
        # only calculate upper triangle of matrix
        for j in range(i + 1, num_of_seq):
            p = calculate_p_distance(sequences[sequences_id[i]],
                                     sequences[sequences_id[j]])
            dist = jukes_cantor_distance(p)
            # print(f"Distance between {sequences_id[i]} and {sequences_id[j]}: p={p:.4f}, d={dist:.4f}")
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    # print("Distance matrix:")
    print(dist_matrix)
    return dist_matrix, sequences_id

import numpy as np

def get_true_positive_rate_under_thresholds(probs, label):
    """
    return a list which included the true posite rate under different threshold, 
    ranges from 0 to 1, step by 0.01.
    """
    probs = np.array(probs)
    label = np.array(label)

    thresholds = np.linspace(0, 1, 101)  # Generate 100 evenly-spaced thresholds between 0 and 1
    positive_probs = []  # Store the true positive rate at each threshold

    for threshold in thresholds:
        # Samples predicted as positive
        predicted_positive_indices = np.where(probs > threshold)[0]

        true_positive_count = np.sum(label[predicted_positive_indices] == 1)

        # Count of samples predicted as positive
        predicted_positive_count = len(predicted_positive_indices)

        # Compute the true positive rate
        true_positive_rate = true_positive_count / predicted_positive_count if predicted_positive_count > 0 else 0
        
        positive_probs.append(true_positive_rate)

    return  positive_probs

"""
how to use
just get_true_positive_rate_under_thresholds( probs, label )
"""


def add_gaussian_noise(matrix, scale_factor=0.35):
    """
    Add Gaussian noise to a matrix. The noise standard deviation is scaled by
    each matrix element multiplied by scale_factor.

    Parameters:
    - matrix: Input matrix (numpy array).
    - scale_factor: Scaling factor controlling noise magnitude. Default is 0.35.

    Returns:
    - noisy_matrix: Matrix with added noise (numpy array).
    """
    noise = np.random.normal(0, scale_factor * np.abs(matrix), matrix.shape)
    noisy_matrix = matrix + noise
    return noisy_matrix

''''
Example Usage:
b = add_gaussian_noise( a )
'''
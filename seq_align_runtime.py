import csv
import matplotlib.pyplot as plt
import time
import random

DEFAULT_GAP_COST = 1  # default gap cost if not defined in the cost file

def read_cost_matrix(filename):
    """
    Reads the cost matrix from the CSV cost file.
    
    Expects a file with a header row. The first cell of the header is a placeholder (e.g., "*")
    and is ignored. The remaining header cells are the symbols.
    
    Each subsequent row should have a row label (a symbol) in the first cell,
    and then cost values (as integers) corresponding to the header symbols.
    
    Returns:
        A dictionary mapping (symbol1, symbol2) -> cost (int)
    """
    cost_dict = {}
    with open(filename, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    
    # The header row: first element is a placeholder, so ignore it.
    header = rows[0]
    symbols = header[1:]
    
    # Process each subsequent row.
    for row in rows[1:]:
        # The first token is the row symbol.
        row_symbol = row[0]
        # The remaining tokens are the cost values corresponding to the header symbols.
        for col_symbol, cost_val in zip(symbols, row[1:]):
            cost_dict[(row_symbol, col_symbol)] = int(cost_val)
    return cost_dict

def get_cost(x, y, cost_dict):
    """
    Retrieves the cost for aligning symbol x with symbol y.
    
    It first checks if (x, y) exists in cost_dict. If not, since the matrix is symmetric,
    it checks for (y, x). Additionally, if one of the symbols is '-' (a gap) and no entry
    is found, it returns a default gap cost.
    
    Args:
        x (str): First symbol.
        y (str): Second symbol.
        cost_dict (dict): Dictionary of alignment costs.
    
    Returns:
        int: The cost for aligning x with y.
    
    Raises:
        KeyError: If no cost is found (and neither symbol is a gap).
    """
    # If either symbol is a gap, try to get the cost; if not defined, return default.
    if x == '-' or y == '-':
        if (x, y) in cost_dict:
            return cost_dict[(x, y)]
        elif (y, x) in cost_dict:
            return cost_dict[(y, x)]
        else:
            return DEFAULT_GAP_COST
    # For non-gap symbols.
    if (x, y) in cost_dict:
        return cost_dict[(x, y)]
    elif (y, x) in cost_dict:
        return cost_dict[(y, x)]
    else:
        raise KeyError(f"Cost for aligning '{x}' and '{y}' not found.")

def align_sequences(seq1, seq2, cost_dict):
    """
    Computes the optimal global alignment for two sequences using dynamic programming.
    
    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        cost_dict (dict): Dictionary with alignment costs for symbol pairs.
    
    Returns:
        A tuple (aligned_seq1, aligned_seq2, total_cost)
    """
    m, n = len(seq1), len(seq2)
    
    # Initialize the DP table: dimensions (m+1) x (n+1)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize the first column (seq1 vs. gaps)
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + get_cost(seq1[i - 1], '-', cost_dict)
    
    # Initialize the first row (gaps vs. seq2)
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + get_cost('-', seq2[j - 1], cost_dict)
    
    # Fill in the DP table.
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost_diag = dp[i - 1][j - 1] + get_cost(seq1[i - 1], seq2[j - 1], cost_dict)
            cost_up   = dp[i - 1][j]     + get_cost(seq1[i - 1], '-', cost_dict)
            cost_left = dp[i][j - 1]     + get_cost('-', seq2[j - 1], cost_dict)
            dp[i][j] = min(cost_diag, cost_up, cost_left)
    
    # Backtrace to recover one optimal alignment.
    aligned_seq1 = []
    aligned_seq2 = []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i - 1][j - 1] + get_cost(seq1[i - 1], seq2[j - 1], cost_dict):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i - 1][j] + get_cost(seq1[i - 1], '-', cost_dict):
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
    
    # Since we built the alignment backwards, reverse the lists.
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    
    alignment_str1 = ''.join(aligned_seq1)
    alignment_str2 = ''.join(aligned_seq2)
    total_cost = dp[m][n]
    
    return alignment_str1, alignment_str2, total_cost

def main(length):
    # Filenames for the cost matrix and input sequences.
    cost_filename = "imp2cost.txt"
    input_filename = "imp2outputruntime.txt"
    output_filename = "imp2outputruntime.txt"

    inputs = []
    for i in range(0,10):
        dnaseq1 = ''.join(random.choices('AGTC', k=length))
        dnaseq2 = ''.join(random.choices('AGTC', k=length))
        inputs.append(f"{dnaseq1}, {dnaseq2}")

    with open(output_filename, 'w') as outfile:
        for result_line in inputs:
            outfile.write(result_line + "\n")
    
    # Read the cost matrix.
    cost_dict = read_cost_matrix(cost_filename)
    
    results = []
    # Process the input file containing sequence pairs.
    with open(input_filename, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue  # Skip empty lines.
            # Each line should contain two sequences separated by a comma.
            parts = line.split(',')
            if len(parts) != 2:
                print(f"Skipping malformed line: {line}")
                continue
            seq1 = parts[0].strip()
            seq2 = parts[1].strip()
            
            # Compute the alignment and total cost.
            aligned_seq1, aligned_seq2, cost_val = align_sequences(seq1, seq2, cost_dict)
            
            # Format the result as: aligned_seq1, aligned_seq2:total_cost
            results.append(f"{aligned_seq1}, {aligned_seq2}:{cost_val}")
    
    # Write the results to the output file.
    with open(output_filename, 'w') as outfile:
        for result_line in results:
            outfile.write(result_line + "\n")

if __name__ == "__main__":
    inputsizes = [500,1000,2000,4000,5000]
    runtimes=[]
    for i in range(0,len(inputsizes)):

        start_time = time.time()
        main(inputsizes[i])
        end_time = time.time()
        runtime = end_time - start_time
        runtimes.append(runtime)

    # Create the plot
    plt.figure(figsize=(10, 6))

    # Plot the empirical runtime
    plt.plot(inputsizes, runtimes, 'o-', color='green', label='Align')

    # Set both axes to logarithmic scale
    """plt.xscale('log')
    plt.yscale('log')"""

    # Labeling the axes
    plt.xlabel('Input Size (n)', fontsize=14)
    plt.ylabel('Runtime (seconds)', fontsize=14)

    # Adding a title
    plt.title('Empirical Runtime to Align', fontsize=16)

    # Adding a legend
    plt.legend(fontsize=12)

    # Adding grid for better readability
    plt.grid(True, which="both", ls="--", linewidth=0.5)

    # Show the plot
    plt.tight_layout()
    plt.show()

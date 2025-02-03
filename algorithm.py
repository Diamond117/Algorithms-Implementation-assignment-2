import numpy as np

def edit_distance_with_backtrace(s, t):
    m, n = len(s), len(t)
    D = np.zeros((m+1, n+1), dtype=int)
    ptr = np.zeros((m+1, n+1), dtype=tuple)
    
    # Initialize base cases
    for i in range(1, m+1):
        D[i][0] = i
        ptr[i][0] = (i-1, 0)
    for j in range(1, n+1):
        D[0][j] = j
        ptr[0][j] = (0, j-1)
    
    # Fill DP table
    for i in range(1, m+1):
        for j in range(1, n+1):
            cost = 0 if s[i-1] == t[j-1] else 1
            choices = [(D[i-1][j] + 1, (i-1, j)),   # Deletion
                       (D[i][j-1] + 1, (i, j-1)),   # Insertion
                       (D[i-1][j-1] + cost, (i-1, j-1))]  # Substitution
            
            D[i][j], ptr[i][j] = min(choices, key=lambda x: x[0])
    #TO DO: write the table as a file

    
    
    # Backtrace to reconstruct alignment
    aligned_s, aligned_t = [], []
    i, j = m, n
    while (i, j) != (0, 0):
        prev_i, prev_j = ptr[i][j]
        if prev_i == i - 1 and prev_j == j - 1:  # Match/Substitution
            aligned_s.append(s[i-1])
            aligned_t.append(t[j-1])
        elif prev_i == i - 1:  # Deletion
            aligned_s.append(s[i-1])
            aligned_t.append("-")
        else:  # Insertion
            aligned_s.append("-")
            aligned_t.append(t[j-1])
        i, j = prev_i, prev_j
    
    return D[m][n], "".join(aligned_s[::-1]), "".join(aligned_t[::-1])

# testing
s1, s2 = "AAATGTGTGTGTTCCCCAACGATGTCTCTAGAAGACGAACATCCC", "ATGGAAACGTGAACCTAACTAACACATATGGATCCGACTGACGTTCTCTGATGTAGCCT"
distance, aligned_s1, aligned_s2 = edit_distance_with_backtrace(s1.replace('-', ''), s2.replace('-', ''))
print(f"Edit Distance: {distance}")
print(f"Alignment:\n{aligned_s1}\n{aligned_s2}")

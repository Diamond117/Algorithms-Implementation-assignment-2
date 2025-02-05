import numpy as np
#import imp2cost # type: ignore

def edit_distance_with_backtrace(s, t):
    m, n = len(s), len(t)
    D = np.zeros((m+1, n+1), dtype=int)
    ptr = np.zeros((m+1, n+1), dtype=tuple)

    #TO DO: read in cost table as referance
    # open the sample file used 
    file = open('imp2cost.txt') 
    
    # read the content of the file opened 
    content = file.readlines() 

    # Initialize base cases
    for i in range(1, m+1):
        columnIndex = content[0].index(s[i-1])
        tablecost = int(content[1][columnIndex])
        D[i][0] = tablecost + D[i-1][0]
        ptr[i][0] = (i-1, 0)
    for j in range(1, n+1):
        columnIndex = content[0].index(t[j-1])
        tablecost = int(content[1][columnIndex])
        D[0][j] = tablecost + D[0][j-1]
        ptr[0][j] = (0, j-1)
    
    
    # Fill DP table
    for i in range(1, m+1):
        for j in range(1, n+1):
            rowIndex = int(content[0].index(s[i-1])/2) 
            columnIndex = content[0].index(t[j-1])
            tablecost = int(content[rowIndex][columnIndex])
            choices = [(D[i-1][j] + tablecost, (i-1, j)),   # left
                       (D[i][j-1] + tablecost, (i, j-1)),   # down
                       (D[i-1][j-1] + tablecost, (i-1, j-1))]  # diagonal
            
            D[i][j], ptr[i][j] = min(choices, key=lambda x: x[0])
    #TO DO: write the table as a file
    for i in range(0, m+1):
        lin=""
        for j in range(0, n+1):
            lin += str(D[m-i][j])
            lin += ", "
        print(lin)
    
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
s1, s2 = "ATGC", "ATGC"
distance, aligned_s1, aligned_s2 = edit_distance_with_backtrace(s1.replace('-', ''), s2.replace('-', ''))
print(f"Edit Distance: {distance}")
print(f"Alignment:\n{aligned_s1}\n{aligned_s2}")

"""*,-,A,T,G,C
    -,0,1,2,1,3
    A,1,0,1,5,1
    T,2,1,0,9,1
    G,1,5,9,0,1
    C,3,1,1,1,0"""
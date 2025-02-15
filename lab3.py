
import numpy as np
from lab2 import MySimplexMethodMainPhase

def MySimplexMethodStartPhase(c: np.ndarray, A: np.ndarray, b: np.ndarray):
# step 1 - vect prav ch b neotr
    negative_index = b < 0
    A[negative_index] *= -1
    b[negative_index] *= -1
    print("\nA")
    print(A)
    print("\nb")
    print(b)

# step 2 - vspomog zad lin progr
    m, n = A.shape
    
    c_wave = np.zeros(n + m)
    c_wave[n:] = -1
    print("\nC_Wave")
    print(c_wave)
    
    # x (n+1..n+m) - iskusstv perem

    # pris sprav ed matr por m
    identity_matrix = np.eye(m)
    A_wave = np.hstack((A, identity_matrix))
    print("\nA_Wave")
    print(A_wave)
    
# step 3 - stroim nachal baz pl
    x_wave = np.zeros(n + m)
    x_wave[n:] = b[:]
    print("\nx_wave")
    print(x_wave)
    
    B = [i + n for i in range(1, m + 1)]
    print("\nB")
    print(B)
    
# step 4 - resh vspm zad
    x_wave, B = MySimplexMethodMainPhase(c_wave, x_wave, A_wave, B)
    print("\nx_wave")
    print(x_wave)
    print("\nB")
    print(B)
    
# step 5 - prov usl sovm usl sovm =0 - sist sovm
    if np.any(x_wave[n] > 0):
        raise Exception("problem is not feasible")
    
# step 6 - formir dop pl x
    x = x_wave[:n]
    print("\nx")
    print(x)
    
    while True:
# step 7 - vozm zaversh (B prinadl 1..n)
        
        if all(bi <= n for bi in B):
            return x, B, A, b
        
# step 8 - viberem  v B max ind issk per
        j_k = max(B)
        k = B.index(j_k) + 1
        print("\nj_k")
        print(j_k)
        
        print("\nk")
        print(k)
        
# step 9 - vich vect dla kagdogo ind
        j_NB = [i+1 for i in range(n) if all(bi - 1 != i for bi in B)]
        print("\nj_NB")
        print(j_NB)
        
        A_B = A_wave[:, np.array(B) - 1]
        A_B_inverted = np.linalg.inv(A_B)
        print("\nA_B_inverted")
        print(A_B_inverted)
        
        l = [(j, A_B_inverted.dot(A_wave[:, j - 1])) for j in j_NB]
        print("\nl")
        print(l)

# step 10 - esl naid ind tak chto ... to zamenim
        found = False
        for j,l_j in l:
            if l_j[k-1] != 0:
                found = True
                B[k - 1] = j

# step 11 esli dla lubogo ind vipoln lk=0
        if not found:
            i = j_k - n
            
            A = np.delete(A, i - 1, axis=0)
            b = np.delete(b, i - 1)
            B = np.delete(B, k - 1)
            A_wave = np.delete(A_wave, i - 1, axis=0)

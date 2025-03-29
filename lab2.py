from hmac import new
from math import inf
import numpy as np
from lab1 import MyMatrixInversion

def MySimplexMethodMainPhase(c: np.ndarray, x: np.ndarray, A: np.ndarray, B: np.ndarray):
    Is_First_Lap = True
    new_index = 0
    new_column = A[new_index]
    m = A.shape[0]
    iteration = 1
    
    while True:
        #STEP 1
        A_B = A[:, np.array(B) - 1]
        
        #print ('A base')
        #print (A_B)
        #print('\n')
        
        if Is_First_Lap:
            A_I_B= np.linalg.inv(A_B)
        else:
            A_I_B= MyMatrixInversion(A_I_B, new_column, new_index)
        
        if (A_I_B is None):
            #print('couldnt inverse matrix')
            return
        
        #print('A base inverted')
        #print(A_I_B)
        #print('\n')
            
        #STEP 2
        c_B = c[np.array(B) - 1]
        
        #print('c component with base indexes')
        #print(c_B)
        #print('\n')
        
        #STEP 3
        u = c_B.dot(A_I_B)
        
        #print('u')
        #print(u)
        #print('\n')
        
        #STEP 4
        delta = u.dot(A) - c
        
        #print('estimates vector')
        #print(delta)
        #print('\n')
        
        #STEP 5
        if np.all(delta >= 0):
            #print('found optimal plan')
            return x, B
        
        #STEP 6
        j0 = (delta < 0).argmax() + 1
        
        #print('index of first negative number')
        #print(j0)
        #print('\n')
        
        #STEP 7
        z = A_I_B.dot(A[:,j0 - 1])
        
        #print('z')
        #print(z)
        #print('\n')
        
        #STEP 8-9(didn't create whole vector because it's unnecessary)
        teta0 = inf
        k = -1
        for i in range(m):
            if z[i] > 0 and teta0 > x[int(B[i] - 1)] / z[i] :
                teta0 = x[int(B[i] - 1)] / z[i]    
                k = i + 1
            else:
                teta0
            
        
        #print('teta0')
        #print(teta0)
        #print('\n')
        
        #STEP 10
        if teta0 == inf:
            #print('objective function is not limited from above on a set of acceptable plans')
            return
        
        #STEP 11
        j_Star = B[k-1]
        
        #print('k')
        #print(k)
        #print('\n')
        #print('j_Star = k')
        #print(j_Star)
        #print('\n')
        
        #STEP 12
        B[k-1] = j0
        
        #STEP 13
        for i in range (m):
            x[int(B[i] - 1)] = x[int(B[i] - 1)] - teta0*z[i] if i != k-1 else x[int(B[i] - 1)]
            
        x[j_Star - 1] = 0
        x[j0 - 1] = teta0
        #print(x)
        #print('\n')
            
        
        
        new_index = k
        new_column = A[:,B[k-1] - 1]
        Is_First_Lap = False
        iteration+=1
    
# Коэффициенты целевой функции
c = np.array([3, 2, 0, 0])

# Матрица ограничений
A = np.array([
    [1, 1, 1, 0],
    [2, 1, 0, 1]
])

# Базисные переменные
B = [3, 4]

# Начальный базисный план
x = np.array([0, 0, 4, 5])

# Запуск симплекс-метода
solution, final_basis = MySimplexMethodMainPhase(c, x, A, B)

# Вывод результатов
print("Оптимальное решение:")
print(solution)
print("\nИтоговый базис:")
print(final_basis)
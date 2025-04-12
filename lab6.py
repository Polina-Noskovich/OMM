import numpy as np
import math

# Задача квадратичного программирования —> это оптимизационная задача -> Допустимый план x -> правильный опорный план
def is_negative(value):
    for item in range(len(value)):
        if value[item] < 0:  # для всяк вектора x выполн нерав
            return value[item], item
    return None, None


def get_submatrix_d(D, J_):
    new_D = [[0] * len(J_) for _ in range(len(J_))]
    i, j = 0, 0
    for first_index in J_:
        for second_index in J_:
            new_D[i][j] = D[first_index][second_index]
            j += 1
        i += 1    # обратимая матрица Аb
        j = 0
    return new_D

# найдется подмножество Jb* множества индексов переменных
# для каждого индекса j ∈ Jbзвездочка выполняется ∆j (x) = 0, где 
# c'(x) = c' + x'*D; u'(x) = -cb'(x)*Ab^-1; ∆'(x) = u'(x)*A + c'(x)

def get_vector_b_star(j0, J_, D, A):
    A_j = A[:, j0]
    D_j0 = D[:, j0]
    D_j = []
    for item in range(len(D_j0)):
        if item in J_:
            D_j.append(D_j0[item])
    b_star = np.concatenate((D_j, A_j), axis=0)
    return b_star

# следующая блочная матрица обратима
# H = (D* -- A(b*)')
#     (A(b*) --   0)
# где D* — это подматрица матрицы D, составленная из элементов,
# стоящих на пересечении строк и столбцов с индексами из множества J(b*);
# A(b*) — матрица, состоящая из столбцов матрицы A с индексами из множества J(b*)

def update(Jb, Jb_star, j0, teta_j0, B, A):
    # случай 1
    if j0 == teta_j0:
        Jb_star.append(teta_j0)
        return
    # случай 2
    Jb_star_not_in_Jb = [i for i in Jb_star if i not in Jb]
    if teta_j0 in Jb_star_not_in_Jb:
        Jb_star.remove(teta_j0)
    # случай 3 и 4
    if teta_j0 in Jb:
        s = Jb.index(teta_j0)
        exists_index = False
        j_pos = -1
        for i in Jb_star_not_in_Jb:
            value = np.dot(B, A[:, i])
            if value.item(s) != 0:
                exists_index = True
                j_pos = i
                break
        if exists_index:
            Jb_star.remove(teta_j0)
            Jb[s] = j_pos
            return
        Jb[s] = j0
        Jb_star[Jb_star.index(teta_j0)] = j0

# Множество Jb называется опорой ограничений, а множество J(b*) — расширенной опорой ограничений.


def catere_vector_l(D, B, Ab, J_, j0):
    l = [0 for _ in range(len(D))]
    l[j0] = 1
    submatrix_D = get_submatrix_d(D, J_)
    matrix_Z = np.zeros((len(J_), len(J_)))
    H = np.bmat([[submatrix_D, B], [Ab, matrix_Z]])
    b_star = get_vector_b_star(j0, J_, D, A)
    x_items = np.array((-1) * np.dot((np.linalg.inv(H)), b_star))[0]
    i = 0
    for item in J_:
        l[item] = x_items[i]
        i += 1
    return l


def get_tetes(beta, j0, delta_x, J_, l):
    tetes = [0 for _ in range(len(D))]
    if beta == 0:
        tetes[j0] = math.inf
    else:
        tetes[j0] = math.fabs(delta_x[j0]) / beta
    for j in J_:
        if l[j] >= 0:
            tetes[j] = math.inf
        else:
            tetes[j] = ((-1) * x[j]) / l[j]
    return tetes


def min_tetes_and_index(tetes):
    teta_min = min(tetes)
    for i in range(len(tetes)):
        if tetes[i] == teta_min:
            return teta_min, i
    return None, None


def quadratic_programming(A, D, b, c, J_on, J_, x):
    while True:
        Ab = A[:, J_on]
        B = np.linalg.inv(Ab)
        c_x = c + np.dot(D, x)
        c_b = np.array([c_x[j] for j in J_on])
        u = (-1) * (np.matmul(c_b, B))
        delta_x = np.matmul(u, A) + c_x
        delta_min, j0 = is_negative(delta_x)
        if j0 is None:
            print('Допустимый план' + str(x))
            print('Опора ограничений' + str(J_on))
            print('Расширенная опора ограничений' + str(J_))
            return x
        l = catere_vector_l(D, B, Ab, J_, j0)
        beta = np.dot((np.dot(l, D)), np.transpose(l))
        tetes = get_tetes(beta, j0, delta_x, J_, l)
        teta_0, teta_j0 = min_tetes_and_index(tetes)
        if teta_0 == math.inf:
            return None
        x = x + np.dot(teta_0, l)
        update(J_on, J_, j0, teta_j0, B, A)


if __name__ == '__main__':
    A = np.array([[2, 1, 0],
                  [1, 0, 1]])
    D = np.array([[6, -2, 0],
                  [-2, 2, 0],
                  [0, 0, 2]])
    b = np.array([2, 1])
    c = np.array([-1, -1, -1])
    J_on = [1, 2]
    J_ast = [1, 2]
    x = np.array([0, 2, 1])
    

    x = quadratic_programming(A, D, b, c, J_on, J_ast, x)
    print(x)


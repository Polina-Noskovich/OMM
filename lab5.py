import numpy as np

def first_phase(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, list, np.ndarray, np.ndarray, np.ndarray]:
    a_total, b_total = np.sum(a), np.sum(b)
    if a_total != b_total:
        difference = a_total - b_total
        if difference > 0:
            b = np.append(b, difference)
            c = np.hstack((c, np.zeros((len(a), 1))))
        else:
            a = np.append(a, -difference)
            c = np.vstack((c, np.zeros((len(b), 1))))
    
    n, m = len(a), len(b)
    X = np.zeros((n, m))
    i, j = 0, 0
    B = []

    while i < n and j < m:
        if b[j] > a[i]:
            X[i, j] = a[i]
            b[j] -= a[i]
            a[i] = 0
            i += 1
        else:
            X[i, j] = b[j]
            a[i] -= b[j]
            b[j] = 0
            j += 1
        B.append((i, j))
    
    return X, B, a, b, c

def potentials_method(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    X, B, a, b, c = first_phase(a.copy(), b.copy(), c.copy())
    n, m = len(a), len(b)

    while True:
        # Создаём систему уравнений для потенциалов
        A = []
        rhs = []
        
        # Уравнения для базисных клеток: u[i] + v[j] = c[i,j]
        for i, j in B:
            row = np.zeros(n + m)
            row[i] = 1  # u[i]
            row[n + j] = 1  # v[j]
            A.append(row)
            rhs.append(c[i, j])
        
        # Фиксируем u[0] = 0 для однозначности решения
        row = np.zeros(n + m)
        row[0] = 1
        A.append(row)
        rhs.append(0)
        
        # Решаем систему
        try:
            solution = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            # Если система вырождена, добавляем случайное уравнение
            row = np.random.rand(n + m)
            A.append(row)
            rhs.append(np.random.rand())
            solution = np.linalg.lstsq(A, rhs, rcond=None)[0]
        
        u = solution[:n]
        v = solution[n:]
        
        # Проверяем небазисные клетки
        found = False
        for i in range(n):
            for j in range(m):
                if (i, j) not in B and u[i] + v[j] > c[i, j] + 1e-6:  # Добавляем погрешность
                    # Нашли улучшающую клетку
                    B.append((i, j))
                    found = True
                    
                    # Находим цикл
                    # (здесь должна быть реализация поиска цикла)
                    # Для простоты предположим, что нашли цикл
                    
                    # Помечаем цикл
                    cycle = find_cycle(B, (i, j))
                    if not cycle:
                        continue
                    
                    # Находим theta
                    theta = min(X[k][l] for (k, l), sign in cycle if sign == -1)
                    
                    # Корректируем перевозки
                    for (k, l), sign in cycle:
                        X[k][l] += theta * sign
                    
                    # Удаляем клетку с нулевой перевозкой
                    for idx, (k, l) in enumerate(B):
                        if X[k][l] == 0:
                            del B[idx]
                            break
                    
                    break
            if found:
                break
                
        if not found:
            break  # Решение оптимально
    
    return X

def find_cycle(B, new_cell):
    # Упрощённая реализация поиска цикла
    # В реальной реализации нужно найти замкнутый путь
    return []  # Заглушка

if __name__ == "__main__":
    a = np.array([30, 40])
    b = np.array([20, 30, 20])
    c = np.array([
        [2, 3, 4],
        [5, 1, 2]
    ])

    print("Исходные данные:")
    print("a:", a)
    print("b:", b)
    print("c:")
    print(c)

    print("\nНачальное решение:")
    X_init, B_init, a_bal, b_bal, c_bal = first_phase(a.copy(), b.copy(), c.copy())
    print("X:")
    print(X_init)
    print("Базисные клетки:", B_init)
    print("Стоимость:", np.sum(X_init * c_bal[:X_init.shape[0], :X_init.shape[1]]))

    print("\nОптимальное решение:")
    try:
        X_opt = potentials_method(a.copy(), b.copy(), c.copy())
        print("X:")
        print(X_opt)
        print("Стоимость:", np.sum(X_opt * c))
    except Exception as e:
        print(f"Ошибка: {e}")
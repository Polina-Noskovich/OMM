import numpy as np

# Первая фаза (метод северо-западного угла)
# Суть метода:
# Мы начинаем заполнять транспортную таблицу с верхнего левого угла (северо-западного).
# В каждой ячейке (i, j) мы заполняем перевозимое количество груза, 
# ориентируясь на доступные запасы у поставщика a[i] и потребности у потребителя b[j].
# Когда запас исчерпан, переходим к следующему поставщику (i+1).
# Когда потребность удовлетворена, переходим к следующему потребителю (j+1).
# Процесс продолжается, пока все грузы не распределены.
def first_phase(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, list, np.ndarray, np.ndarray, np.ndarray]:
    # Балансировка запасов и потребностей
    a_total, b_total = np.sum(a), np.sum(b)  # Считаем общий запас и потребность
    if a_total != b_total:
        difference = a_total - b_total
        if difference > 0:
            # Добавляем фиктивного потребителя с нулевой стоимостью
            b = np.append(b, difference)
            c = np.hstack((c, np.zeros((len(b), 1))))
        else:
            # Добавляем фиктивного поставщика с нулевой стоимостью
            a = np.append(a, -difference)
            c = np.vstack((c, np.zeros((len(a), 1))))
    
    # Определение размеров матрицы
    n, m = len(a), len(b)
    X = np.zeros((n, m))  # Инициализация таблицы перевозок нулями
    i, j = 0, 0  # Индексы текущего поставщика и потребителя
    B = []  # Список базисных клеток

    # Заполнение методом северо-западного угла
    while i < n and j < m:
        if b[j] > a[i]:
            X[i, j] = a[i]  # Отправляем весь запас поставщика
            b[j] -= a[i]  # Уменьшаем потребность
            a[i] = 0  # Запас поставщика исчерпан
            i += 1  # Переход к следующему поставщику
        else:
            X[i, j] = b[j]  # Отправляем весь объем, который требуется потребителю
            a[i] -= b[j]  # Уменьшаем запас поставщика
            b[j] = 0  # Потребность удовлетворена
            j += 1  # Переход к следующему потребителю
        B.append((i, j))  # Запоминаем базисные клетки
    
    return X, B, a, b, c

# Метод потенциалов
# Суть метода:
# После нахождения базисного решения (методом северо-западного угла)
# оцениваем оптимальность решения с помощью потенциалов.
# Для этого:
# 1. Вводятся потенциалы u[i] (для поставщиков) и v[j] (для потребителей).
# 2. Решается система уравнений u[i] + v[j] = c[i, j] для всех базисных клеток.
# 3. Проверяются оценки для всех небазисных клеток: если оценка >= стоимости c[i, j], решение оптимально.
# 4. Если найдена клетка с отрицательной оценкой, строится замкнутый цикл, 
#    выполняется пересчет перевозок, и процесс повторяется.
def potentials_method(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    # Получение начального базисного решения
    X, B, a, b, c = first_phase(a, b, c)
    n, m = len(a), len(b)

    while True:
        x, y = [], []  # Матрицы для системы уравнений потенциалов
        for i, j in B:
            u, v = [0] * n, [0] * m
            u[i], v[j] = 1, 1  # Формируем уравнение u[i] + v[j] = c[i,j]
            x.append(u + v)
            y.append(c[i, j])

        x.append([1] + [0] * (n + m - 1))  # Добавляем уравнение для фиксирования одного потенциала
        y.append(0)

        # Решаем систему уравнений методом линейной алгебры
        result = np.linalg.solve(x, y)
        u, v = result[:n], result[n:]

        new_position = None
        for i in range(n):
            for j in range(m):
                if (i, j) not in B and u[i] + v[j] > c[i, j]:
                    new_position = (i, j)  # Ищем выгодную клетку
                    break
            if new_position:
                break

        if not new_position:
            return X  # Если нет выгодных клеток, решение оптимально
        
        B.append(new_position)  # Добавляем новую клетку в базис

        # Находим базисный цикл
        corner_B = B[:]
        while True:
            changes = False
            for i in range(n):
                connected = [j for j in range(m) if (i, j) in corner_B]
                if len(connected) <= 1:
                    for j in connected:
                        corner_B.remove((i, j))
                        changes = True

            for j in range(m):
                connected = [i for i in range(n) if (i, j) in corner_B]
                if len(connected) <= 1:
                    for i in connected:
                        corner_B.remove((i, j))
                        changes = True

            if not changes:
                break

        # Помечаем клетки знаками + и -
        marked_B = {position: None for position in corner_B}
        marked_B[new_position] = True
        add_plus_or_minus(new_position, marked_B)

        # Определяем минимальное значение для корректировки перевозок
        theta = min(X[i, j] for i, j in marked_B if not marked_B[(i, j)])
        for i, j in marked_B:
            if marked_B[(i, j)]:
                X[i, j] += theta
            else:
                X[i, j] -= theta

        # Удаляем клетку с минимальным значением из базиса
        B.remove(next((pos for pos in marked_B if not marked_B[pos]), None))

# Функция для маркировки базисного цикла знаками + и -
def add_plus_or_minus(position: tuple[int, int], B: dict) -> None:
    for i, j in list(B.keys()):
        if (position[0] == i or position[1] == j) and B[(i, j)] is None:
            B[(i, j)] = not B[position]  # Чередуем знаки + и -
            add_plus_or_minus((i, j), B)

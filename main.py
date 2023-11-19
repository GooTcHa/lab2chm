import math
from copy import deepcopy


def fill_matrix(A: list, b: list, x: list, r: list, n: int):
    #Заполняем матрицы A и b
    s = 0.0
    s_line = 0.0
    diagonal_dominance = True
    for i in range(n):
        x.append(9 + i)
    for i in range(n):
        #Проверяем диагональное приобладание матрицы А
        if not diagonal_dominance:
            raise Exception('Matrix doesn\'t have diagonal domination')
        A.append([])
        for j in range(n):
            if i == j:
                A[i].append(11 * math.sqrt(i + 1))
            else:
                A[i].append(0.001 * (i + 1) / (j + 1))
                #Сумма элементов не лежащих на главной диагонали
                s_line += math.fabs(A[i][j])
            s += A[i][j] * x[j]
        r.append(s_line)
        #Проверям усовие строгово диагонального приабладания
        if s_line > A[i][i]:
            diagonal_dominance = False
        s_line = 0.0
        b.append(s)
        s = 0


def solve_equations_system_by_Jacobi(A: list, b: list, x: list, eps: float, k: int,  n: int):
    x0 = deepcopy(b)
    x1 = deepcopy(b)
    comparison_vector: list
    answer = []
    final_k = 0
    #итерации
    for l in range(k):
        for i in range(n):
            for j in range(n):
                if i != j:
                    #Подсчитываем X_i^(k+1)
                    x1[i] -= A[i][j] * x0[j]
            # Подсчитываем X_i^(k+1)
            x1[i] /= A[i][i]
        comparison_vector = deepcopy(x1)
        #Находим ||x - x^*||_inf
        for i in range(n):
            comparison_vector[i] -= x0[i]
            comparison_vector[i] = math.fabs(comparison_vector[i])
        #Запоминаем ответ
        answer = x1
        #Запоминаем номер итерации
        final_k = l + 1
        #Проверяем ||x - x^*||_inf < eps
        if max(comparison_vector) < eps:
            break
        x0 = deepcopy(x1)
        x1 = deepcopy(b)
    check_arr= deepcopy(answer)
    #Проверяем погрешность ответа
    for i in range(n):
        check_arr[i] -= x[i]
        check_arr[i] = math.fabs(check_arr[i])
    #Проверяем не достигли ли мы максимального k
    if final_k == k:
        print('k достигло своего максимума и итерации завершились!')
    print(f'На k = {final_k} Решение методом Якоби:\n{answer}\nПроверка решения ||x - x^*||_inf/||x^*||_inf = {max(check_arr)/max(x)}')


def solve_equations_system_by_Gauss_Zeidel(A: list, b: list, x: list, eps: float, k: int,  n: int):
    #Всё аналогично методу якоби за исключением 2 строк
    x0 = deepcopy(b)
    x1 = deepcopy(b)
    comparison_vector: list
    answer = []
    final_k = 0
    for l in range(k):
        for i in range(n):
            for j in range(n):
                if i != j:
                    #Если j меньше i, то x_i^k мы берём как x_i^(k+1)
                    #В противном случае, оставляем x_i^k
                    if i > j:
                        x1[i] -= A[i][j] * x1[j]
                    else:
                        x1[i] -= A[i][j] * x0[j]
            x1[i] /= A[i][i]
        comparison_vector = deepcopy(x1)
        for i in range(n):
            comparison_vector[i] -= x0[i]
            comparison_vector[i] = math.fabs(comparison_vector[i])
        answer = x1
        final_k = l + 1
        if max(comparison_vector) < eps:
            break
        x0 = deepcopy(x1)
        x1 = deepcopy(b)
    check_arr = deepcopy(answer)
    for i in range(n):
        check_arr[i] -= x[i]
        check_arr[i] = math.fabs(check_arr[i])
    if final_k == k:
        print('k достигло своего максимума и итерации завершились!')
    print(f'\nНа k = {final_k} решение методом Гауссa-Зейделя:\n{answer}\nПроверка решения ||x - x^*||_inf/||x^*||_inf = '
          f'{max(check_arr) / max(x)}')


def solve_equations_system_by_Relaxation(A: list, b: list, x: list, eps: float, k: int, w: float,  n: int):
    #Всё аналогично методу Гаусса-Зейделя за исключением w
    x0 = deepcopy(b)
    x1 = deepcopy(b)
    comparison_vector: list
    answer = []
    final_k = 0
    for l in range(k):
        for i in range(n):
            for j in range(n):
                if i != j:
                    if i > j:
                        x1[i] -= A[i][j] * x1[j]
                    else:
                        x1[i] -= A[i][j] * x0[j]
            #Решение полученное методом гаусса зейделя домножаем на w
            x1[i] *= w
            x1[i] /= A[i][i]
            #прибавляем (1-w)x_i^k
            x1[i] += (1 - w)*x0[i]
        comparison_vector = deepcopy(x1)
        for i in range(n):
            comparison_vector[i] -= x0[i]
            comparison_vector[i] = math.fabs(comparison_vector[i])
        answer = x1
        final_k = l + 1
        if max(comparison_vector) < eps:
            break
        x0 = deepcopy(x1)
        x1 = deepcopy(b)
    check_arr = deepcopy(answer)
    for i in range(n):
        check_arr[i] -= x[i]
        check_arr[i] = math.fabs(check_arr[i])
    if final_k == k:
        print('k достигло своего максимума и итерации завершились!')
    print(f'\nНа k = {final_k} решение методом Релаксации c w = {w}:\n{answer}\nПроверка решения ||x - x^*||_inf/||x^*||_inf = '
          f'{max(check_arr) / max(x)}')


def change_A(A: list, b: list, x: list, r: list, n: int):
    for i in range(n):
        A[i][i] = r[i]
        for j in range(n):
            b[i] += A[i][j] * x[j]


if __name__ == '__main__':
    A = []
    b = []
    x = []
    r = []
    n = 5
    #Заполняю матрицу A, b, x, r, где x - точное решение, r - вектор, который будет использоваться для дальнейшего
    #изменения матрицы А
    fill_matrix(A, b, x, r, n)
    print('Матрица А:', *A, '\nМатрица b:', b, '\n', sep='\n')
    #Задаём параметры эпсилон и k
    eps = 0.001
    k = 1000
    #Решаем систему различными способами
    solve_equations_system_by_Jacobi(A, b, x, eps, k, n)
    solve_equations_system_by_Gauss_Zeidel(A, b, x, eps, k, n)
    solve_equations_system_by_Relaxation(A, b, x, eps, k, 0.5, n)
    solve_equations_system_by_Relaxation(A, b, x, eps, k, 1.5, n)

    b = [0] * n
    change_A(A, b, x, r, n)
    print('\nМатрица А после изменения: ', *A, '\nМатрица b после изменения:', b, '\nРешения с изменённой матрицей А:', sep='\n')

    solve_equations_system_by_Jacobi(A, b, x, eps, k, n)
    solve_equations_system_by_Gauss_Zeidel(A, b, x, eps, k, n)
    solve_equations_system_by_Relaxation(A, b, x, eps, k, 0.5, n)
    solve_equations_system_by_Relaxation(A, b, x, eps, k, 1.5, n)
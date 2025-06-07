import math
import matplotlib.pyplot as plt
import numpy as np

def linear_regression(x, y):
    """Выполняет линейную регрессию y = a*x + b и возвращает a, b"""
    N = len(x)
    sum_x = sum(x)
    sum_y = sum(y)
    sum_x2 = sum(xi ** 2 for xi in x)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))

    denominator = N * sum_x2 - sum_x ** 2
    a = (N * sum_xy - sum_x * sum_y) / denominator
    b = (sum_y * sum_x2 - sum_x * sum_xy) / denominator

    return a, b

def get_scu(a, b, x, y):
    """Вычисляет сумму квадратов отклонений (СКУ)"""
    scu_elements = [(1 / (a * xi + b) - yi) ** 2 for xi, yi in zip(x, y)]
    return sum(scu_elements)

def get_sco(scu, N):
    """Вычисляет среднюю ошибку аппроксимации (СКО)"""
    return math.sqrt(scu / N)

def plot_result(x, y, a, b):
    """Строит график исходных данных и аппроксимации"""
    x_plot = np.linspace(min(x) * 0.9, max(x) * 1.1, 300)
    y_plot = 1 / (a * x_plot + b)

    plt.scatter(x, y, color='blue', label='Исходные точки')
    plt.plot(x_plot, y_plot, color='red', label=f'Аппроксимация: $y = 1/({a:.3f}x + {b:.3f})$')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Аппроксимация функцией $y = 1/(a x + b)$')
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    # Исходные данные
    x = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    y = [2.07, 2.41, 1.77, 1.98, 2.24, 2.73, 2.91, 2.56, 3.07]

    # Линеаризация: преобразуем y -> 1/y
    y_transformed = [1 / yi for yi in y]

    # Находим коэффициенты a и b для 1/y = a*x + b
    a, b = linear_regression(x, y_transformed)

    print(f"Коэффициенты: a = {a:.6f}, b = {b:.6f}")
    print(f"Итоговая функция: y = 1 / ({a:.3f} * x + {b:.3f})")

    scu = get_scu(a, b, x, y)
    sco = get_sco(scu, len(x))
    print(f"Средняя квадратическая ошибка (СКО): {sco:.6f}")
    plot_result(x, y, a, b)

if __name__ == "__main__":
    main()
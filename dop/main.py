import cmath
import matplotlib.pyplot as plt
import numpy as np


def f(x):
    return np.cos(2 * np.pi * x)


def generate_sample_points(interval, num_points):
    x_values = np.arange(num_points) * interval / num_points
    y_values = f(x_values)
    return x_values, y_values


def FastFurye(data):
    N = len(data)
    if N <= 1:
        return data
    even_values = FastFurye(data[::2])
    odd_values = FastFurye(data[1::2])
    T = [cmath.exp(-2j * cmath.pi * k / N) * odd_values[k] for k in range(N//2)]
    return [even_values[k] + T[k] for k in range(N//2)] + [even_values[k] - T[k] for k in range(N//2)]


def compute_fourier_coefficients(y_values, num_points):
    coefficients = FastFurye(y_values)
    return np.array(coefficients) / num_points


def trig_interpolation(coefficients, x_eval, interval):
    total = 0.0
    N = len(coefficients)
    for k in range(N):
        freq = k if k <= N//2 else k - N
        total += coefficients[k] * cmath.exp(2j * cmath.pi * freq * x_eval / interval)
    return total.real


def evaluate_interpolation(coefficients, interval, num_points, N):
    x_eval_values = (np.arange(N) + 0.5) * interval / num_points
    y_interp_values = np.array([trig_interpolation(coefficients, x, interval) for x in x_eval_values])
    y_eval_values = f(x_eval_values)
    return x_eval_values, y_interp_values, y_eval_values


def calculate_errors(y_eval, y_interp, x_eval_values):
    abs_errors = np.abs(y_interp - y_eval)
    squared_errors = abs_errors**2
    sum_sq_errors = np.sum(squared_errors)
    rmse = np.sqrt(sum_sq_errors / len(squared_errors))
    max_error = np.max(abs_errors)
    max_error_index = np.argmax(abs_errors)
    x_at_max_error = x_eval_values[max_error_index]
    y_original_at_max = y_eval[max_error_index]
    y_interp_at_max = y_interp[max_error_index]
    return abs_errors, rmse, max_error, x_at_max_error, y_original_at_max, y_interp_at_max


def print_results(y_interp, y_eval, errors):
    abs_errors, rmse, max_error, x_at_max_error, y_orig, y_interp_at_max = errors
    print('Значения интерполированной функции в средних точках разбиений:')
    print(y_interp)
    print('\nЗначения оригинальной функции в средних точках разбиений:')
    print(y_eval)
    print('\nЗначения погрешностей в средних точках разбиений:')
    print(abs_errors)
    print(f'\nМаксимальная погрешность: {max_error:.3e}')
    print(f'При x = {x_at_max_error:.5f}:')
    print(f'  Оригинальная функция: {y_orig:.5f}')
    print(f'  Интерполированная функция: {y_interp_at_max:.5f}')
    print(f'Среднеквадратичная ошибка: {rmse:.3e}')


def plot_results(x_values, y_original, y_interpolated):
    plt.figure(figsize=(14, 7))
    plt.plot(x_values, y_original, 'b-', label='Оригинальная функция', linewidth=1.5)
    plt.plot(x_values, y_interpolated, 'r-', label='Интерполированная функция', linewidth=1.2, alpha=0.7)
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.legend()
    plt.grid(True)
    plt.show()


def main():
    interval = 2
    num_points = 256
    N = 256

    x_values, y_values = generate_sample_points(interval, num_points)
    coefficients = compute_fourier_coefficients(y_values, num_points)
    x_eval, y_interp, y_eval = evaluate_interpolation(coefficients, interval, num_points, N)
    errors = calculate_errors(y_eval, y_interp, x_eval)
    print_results(y_interp, y_eval, errors)
    plot_results(x_eval, y_eval, y_interp)


if __name__ == "__main__":
    main()
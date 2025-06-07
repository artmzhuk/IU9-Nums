import cmath
import matplotlib.pyplot as plt
import numpy as np


def f(x):
    return np.sign(np.cos(2 * np.pi * x))


def generate_sample_points(interval, num_points):
    x_values = [i * interval / num_points for i in range(num_points)]
    y_values = [f(x) for x in x_values]
    return x_values, y_values


def FastFurye(data):
    N = len(data)
    if N <= 1:
        return data
    even_values = FastFurye(data[0::2])
    odd_values = FastFurye(data[1::2])
    factor = cmath.exp(-2j * cmath.pi / N)
    return [even_values[k] + factor ** k * odd_values[k] for k in range(N // 2)] + [
        even_values[k] - factor ** k * odd_values[k] for k in range(N // 2)]


def compute_fourier_coefficients(y_values, num_points):
    coefficients = FastFurye(y_values)
    return [c / num_points for c in coefficients]


def trig_interpolation(coefficients, x_eval, interval):
    result = sum(c * cmath.exp(2j * cmath.pi * n * x_eval / interval)
                 for n, c in enumerate(coefficients))
    return result.real


def evaluate_interpolation(coefficients, interval, num_points, N):
    x_eval_values = [(i * interval + 0.5) / num_points for i in range(N)]
    y_interp_values = [trig_interpolation(coefficients, x, interval) for x in x_eval_values]
    y_eval_values = [f(x) for x in x_eval_values]
    return x_eval_values, y_interp_values, y_eval_values


def calculate_errors(y_eval, y_interp):
    return [y_e - y_i for y_e, y_i in zip(y_eval, y_interp)]


def print_results(y_interp, y_eval, errors):
    print('Значения интерполированной функции в средних точках разбиений:')
    print(y_interp)
    print('\nЗначения оригинальной функции в средних точках разбиений:')
    print(y_eval)
    print('\nЗначения погрешностей в средних точках разбиений:')
    print(errors)
    print(f'\nМаксимальное значение погрешности: {max(errors)}')


def plot_results(x_values, y_original, y_interpolated):
    plt.figure(figsize=(10, 5))
    plt.plot(x_values, y_original, label='Оригинальная функция')
    plt.plot(x_values, y_interpolated, label='Интерполируемая функция')
    plt.legend()
    plt.show()


def main():
    interval = 10
    num_points = 128
    N = 128

    x_values, y_values = generate_sample_points(interval, num_points)
    coefficients = compute_fourier_coefficients(y_values, num_points)
    x_eval, y_interp, y_eval = evaluate_interpolation(coefficients, interval, num_points, N)
    errors = calculate_errors(y_eval, y_interp)
    print_results(y_interp, y_eval, errors)
    plot_results(x_eval, y_eval, y_interp)


if __name__ == "__main__":
    main()
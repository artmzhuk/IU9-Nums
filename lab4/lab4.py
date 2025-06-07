import numpy as np
import matplotlib.pyplot as plt


class RootFinder:
    def __init__(self, tolerance=0.001):
        self.tolerance = tolerance

    def evaluate_function(self, x_value):
        return 2 * x_value ** 3 + 9 * x_value ** 2 - 10

    def evaluate_derivative(self, x_value):
        return 6 * x_value ** 2 + 18 * x_value

    def find_root_bisection(self, left_bound, right_bound):
        if self.evaluate_function(left_bound) * self.evaluate_function(right_bound) >= 0:
            return None, 0, None
        iterations = 0
        while (right_bound - left_bound) / 2 >= self.tolerance:
            iterations += 1
            midpoint = (left_bound + right_bound) / 2
            if self.evaluate_function(left_bound) * self.evaluate_function(midpoint) < 0:
                right_bound = midpoint
            else:
                left_bound = midpoint
        final_root = (left_bound + right_bound) / 2
        actual_error = abs(right_bound - left_bound) / 2
        return final_root, iterations, actual_error

    def find_root_newton(self, initial_guess, max_iterations=100):
        current_x = initial_guess
        previous_x = initial_guess + 2 * self.tolerance  # Инициализация для первой итерации

        for iteration in range(max_iterations):
            func_value = self.evaluate_function(current_x)
            deriv_value = self.evaluate_derivative(current_x)

            # условие 1: Проверка смены знака в ε-окрестности
            test_x = current_x + np.sign(current_x - previous_x) * self.tolerance
            if func_value * self.evaluate_function(test_x) < 0:
                actual_error = abs(func_value)
                return current_x, iteration + 1, actual_error

            if abs(func_value) < self.tolerance:
                actual_error = abs(func_value)
                return current_x, iteration + 1, actual_error

            if deriv_value == 0:
                return None, iteration + 1, None

            previous_x = current_x
            current_x = current_x - func_value / deriv_value

        return None, max_iterations, None



    def locate_all_roots(self, lower_bound, upper_bound, search_step=0.5):
        found_roots = []
        iterations_info = []
        errors_info = []
        current_x = lower_bound
        while current_x <= upper_bound:
            if (self.evaluate_function(current_x) *
                    self.evaluate_function(current_x + search_step) <= 0):
                root, iterations, error = self.find_root_bisection(current_x, current_x + search_step)
                if root is not None and not any(abs(root - r) < self.tolerance for r in found_roots):
                    found_roots.append(root)
                    iterations_info.append(iterations)
                    errors_info.append(error)
            current_x += search_step
        return found_roots, iterations_info, errors_info

    def visualize_results(self, lower_bound, upper_bound, roots_list):
        x_values = np.linspace(lower_bound, upper_bound, 400)
        y_values = self.evaluate_function(x_values)

        plt.figure(figsize=(10, 6))
        plt.plot(x_values, y_values, label="f(x) = 2x³ + 9x² + 10")
        plt.axhline(0, color="black", linewidth=0.5)

        plt.scatter(roots_list, [0] * len(roots_list),
                    color="red", label="Найденные корни", zorder=5)

        plt.title(f"График функции с отмеченными корнями (требуемая точность ε = {self.tolerance})")
        plt.xlabel("x")
        plt.ylabel("f(x)")
        plt.ylim(-20, 20)
        plt.xlim(-20, 20)
        plt.legend()
        plt.grid()
        plt.show()


def execute_root_finding():
    tolerance = 0.001
    finder = RootFinder(tolerance)
    range_start = -10
    range_end = 10

    print(f"Требуемая точность ε = {tolerance}\n")

    bisection_roots, bisection_iterations, bisection_errors = finder.locate_all_roots(range_start, range_end)
    print("Результаты поиска корней методом деления пополам:")
    for idx, (root, iters, error) in enumerate(zip(sorted(bisection_roots), bisection_iterations, bisection_errors), 1):
        print(
            f"Корень #{idx}: x = {root:.5f}, f(x) = {finder.evaluate_function(root):.5f}, итераций: {iters}, достигнутая погрешность: {error:.5f}")

    newton_roots = []
    newton_iterations = []
    newton_errors = []
    for root in bisection_roots:
        newton_root, iters, error = finder.find_root_newton(root)
        if newton_root is not None:
            newton_roots.append(newton_root)
            newton_iterations.append(iters)
            newton_errors.append(error)

    print("\n методом Ньютона:")
    for idx, (root, iters, error) in enumerate(zip(sorted(newton_roots), newton_iterations, newton_errors), 1):
        print(
            f"Корень #{idx}: x = {root:.5f}, f(x) = {finder.evaluate_function(root):.5f}, итераций: {iters}, достигнутая погрешность: {error:.5f}")

    finder.visualize_results(range_start, range_end, bisection_roots)


if __name__ == "__main__":
    execute_root_finding()
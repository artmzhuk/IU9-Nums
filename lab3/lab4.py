import numpy as np
import matplotlib.pyplot as plt


class RootFinder:
    def __init__(self, tolerance=0.001):
        self.tolerance = tolerance

    def evaluate_function(self, x_value):
        return 2 * x_value**3 + 9 * x_value**2 + 10

    def evaluate_derivative(self, x_value):
        return 6 * x_value**2 + 18 * x_value

    def find_root_bisection(self, left_bound, right_bound):
        if self.evaluate_function(left_bound) * self.evaluate_function(right_bound) >= 0:
            return None

        while (right_bound - left_bound) >= self.tolerance:
            midpoint = (left_bound + right_bound) / 2
            if abs(self.evaluate_function(midpoint)) < self.tolerance:
                return midpoint
            if self.evaluate_function(left_bound) * self.evaluate_function(midpoint) < 0:
                right_bound = midpoint
            else:
                left_bound = midpoint
        return (left_bound + right_bound) / 2

    def find_root_newton(self, initial_guess, max_iterations=100):
        current_x = initial_guess
        for _ in range(max_iterations):
            func_value = self.evaluate_function(current_x)
            if abs(func_value) < self.tolerance:
                return current_x
            deriv_value = self.evaluate_derivative(current_x)
            if deriv_value == 0:
                return None
            current_x = current_x - func_value / deriv_value
        return current_x

    def locate_all_roots(self, lower_bound, upper_bound, search_step=0.5):
        found_roots = []
        current_x = lower_bound
        while current_x <= upper_bound:
            if (self.evaluate_function(current_x) * 
                self.evaluate_function(current_x + search_step) <= 0):
                root = self.find_root_bisection(current_x, current_x + search_step)
                if root is not None and not any(abs(root - r) < self.tolerance for r in found_roots):
                    found_roots.append(root)
            current_x += search_step
        return found_roots

    def visualize_results(self, lower_bound, upper_bound, roots_list):
        x_values = np.linspace(lower_bound, upper_bound, 400)
        y_values = self.evaluate_function(x_values)

        plt.figure(figsize=(10, 6))
        plt.plot(x_values, y_values, label="f(x) = 2x³ + 9x² + 10")
        plt.axhline(0, color="black", linewidth=0.5)

        plt.scatter(roots_list, [0] * len(roots_list), 
                   color="red", label="Найденные корни", zorder=5)

        plt.title("График функции с отмеченными корнями")
        plt.xlabel("x")
        plt.ylabel("f(x)")
        plt.legend()
        plt.grid()
        plt.show()


def execute_root_finding():
    finder = RootFinder()
    range_start = -10
    range_end = 10

    bisection_roots = finder.locate_all_roots(range_start, range_end)
    print("Результаты поиска корней методом деления пополам:")
    for idx, root in enumerate(sorted(bisection_roots), 1):
        print(f"Корень #{idx}: x = {root:.5f}, f(x) = {finder.evaluate_function(root):.5f}")

    newton_roots = [finder.find_root_newton(root) for root in bisection_roots]
    print("\nУточнённые значения корней методом Ньютона:")
    for idx, root in enumerate(sorted(newton_roots), 1):
        print(f"Корень #{idx}: x = {root:.5f}, f(x) = {finder.evaluate_function(root):.5f}")

    finder.visualize_results(range_start, range_end, bisection_roots)


if __name__ == "__main__":
    execute_root_finding()
import numpy as np

p = 0.0
q = 4.0
y0 = 3.0
y_prime0 = 4.0
E = 0.001


def f(x):
    return 8.0


def system_of_equations(x, z):
    y, y_prime = z
    dy_dx = y_prime
    dy_prime_dx = -p * y_prime - q * y + f(x)
    return np.array([dy_dx, dy_prime_dx])


def runge_kutta_4_step(f, x, z, h):
    k1 = f(x, z)
    k2 = f(x + h / 2, z + h / 2 * k1)
    k3 = f(x + h / 2, z + h / 2 * k2)
    k4 = f(x + h, z + h * k3)
    return z + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def exact_solution(x):
    return np.cos(2 * x) + 2 * np.sin(2 * x) + 2


def solve_numerically(p, q, y0, y_prime0, E, x_end=1.0):
    h = 0.1  # Фиксированный шаг
    x = 0.0
    z = np.array([y0, y_prime0])
    x_values = [x]
    y_values = [z[0]]
    exact_values = [exact_solution(x)]

    while x < x_end:
        if x + h > x_end:
            h = x_end - x

        z = runge_kutta_4_step(system_of_equations, x, z, h)  # просто один шаг
        x += h

        x_values.append(x)
        y_values.append(z[0])
        exact_values.append(exact_solution(x))

    return x_values, y_values, exact_values



x_values, y_values, exact_values = solve_numerically(p, q, y0, y_prime0, E)

print(" x       | Численное y | Точное y   | Погрешность")
print("---------|-------------|------------|------------")
for x, y_num, y_exact in zip(x_values[:-1], y_values[:-1], exact_values[:-1]):
    error = abs(y_num - y_exact)
    print(f"{x:.5f} | {y_num:.7f} | {y_exact:.7f} | {error:}")
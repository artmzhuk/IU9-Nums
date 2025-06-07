import math

eps = 0.001

def f(x1, x2):
    return 2*x1**2 + 4*x1*x2 + 4*x2**2 - 8*x2 + 0.1*math.log(1 + x1**2*x2**2)

def f_x1(x1, x2):
    return 4*x1 + 4*x2 + 0.1*(2*x1*x2**2)/(1 + x1**2*x2**2)

def f_x2(x1, x2):
    return 4*x1 + 8*x2 - 8 + 0.1*(2*x1**2*x2)/(1 + x1**2*x2**2)

def f_x1_x1(x1, x2):
    denom = (1 + x1**2*x2**2)**2
    term = (2*x2**2*(1 + x1**2*x2**2) - 2*x1**2*x2**4)/denom
    return 4 + 0.1*term

def f_x1_x2(x1, x2):
    denom = (1 + x1**2*x2**2)**2
    term = (4*x1*x2*(1 + x1**2*x2**2) - 4*x1**3*x2**3)/denom
    return 4 + 0.1*term

def f_x2_x2(x1, x2):
    denom = (1 + x1**2*x2**2)**2
    term = (2*x1**2*(1 + x1**2*x2**2) - 2*x1**4*x2**2)/denom
    return 8 + 0.1*term

k = 0
xk, yk = 0.0, 0.0  # Начальное приближение

while max(abs(f_x1(xk, yk)), abs(f_x2(xk, yk))) >= eps:
    phi1 = -(f_x1(xk, yk))**2 - (f_x2(xk, yk))**2
    phi2 = (f_x1_x1(xk, yk) * (f_x1(xk, yk))**2 +
            2 * f_x1_x2(xk, yk) * f_x1(xk, yk) * f_x2(xk, yk) +
            f_x2_x2(xk, yk) * (f_x2(xk, yk))**2)
    t_star = -phi1 / phi2
    xk = xk - t_star * f_x1(xk, yk)
    yk = yk - t_star * f_x2(xk, yk)
    k += 1

print(f'Количество итераций: {k}')
print(f'Минимум: ({xk}, {yk})')
print(f'Значение функции в минимуме: {f(xk, yk)}')
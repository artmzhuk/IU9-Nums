import math
from scipy.integrate import quad

EPS = 0.001

def func(x):
    return 4 * math.atan(x) / (x * x + 1)
 

def trap_method(f, a, b):
    n_arr, I_h_arr, R_arr = [], [], []
    n, R, I_h = 1, 1, 0
    while not (abs(R) < EPS):
        n *= 2
        I_h2 = I_h
        h = (b - a) / n
        s = sum(f(a + i * h) for i in range(1, n))
        I_h = h * ((f(a) + f(b)) / 2 + s)
        R = (I_h - I_h2) / 3 if I_h2 != 0 else 1
        n_arr.append(n)
        I_h_arr.append(I_h)
        R_arr.append(R)
    return n_arr, I_h_arr, R_arr

def simp_method(f, a, b):
    n_arr, I_h_arr, R_arr = [], [], []
    n, R, I_h = 1, 1, 0
    while not (abs(R) < EPS):
        n *= 2
        I_h2 = I_h
        h = (b - a) / n
        s1 = sum(f(a + i * h) for i in range(1, n + 1))
        s2 = sum(f(a + (i - 0.5) * h) for i in range(1, n + 1))
        s3 = sum(f(a + (i - 1) * h) for i in range(1, n + 2))
        I_h = h / 6 * (s1 + 4 * s2 + s3)
        R = (I_h - I_h2) / 15 if I_h2 != 0 else 1
        n_arr.append(n)
        I_h_arr.append(I_h)
        R_arr.append(R)
    return n_arr, I_h_arr, R_arr

if __name__ == "__main__":
    n_t, I_t, R_t = trap_method(func, 0, 4)
    n_s, I_s, R_s = simp_method(func, 0, 4)
    I, err = quad(func, 0, 4)

    print("EPS=",EPS, "\t\tI=",I)

    print("1Метод: \tТрапеции \t\tСимпсона")
    print("n \t\t", n_t[-1],"\t\t\t\t", n_s[-1])
    print("I*\t\t", I_t[-1], "\t\t", I_s[-1])
    print("R \t\t", R_t[-1], "\t", R_s[-1])
    print("I*+R\t\t", I_t[-1] + R_t[-1], "\t\t", I_s[-1] + R_s[-1])

from time import process_time


def function_1(x):
    return x ** 2 - 5


def function_2(x):
    return x ** 2 - 2 * x + 1


def derivative(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2 * h)


def same_sign(a, b):
    return a * b > 0


def bisect(f, x1, x2, eps):
    assert not same_sign(x1, x2)
    mid = x1
    it = 0
    while abs(x1 - x2) > eps:
        if abs(f(mid)) <= 1e-8:
            return mid, it
        elif not same_sign(f(x1), f(mid)):
            x2 = mid
        else:
            x1 = mid
        mid = (x1 + x2) / 2.0
        it += 1
    return mid, it


def secant(f, x0, x1, eps, n_max=1000):
    it = 0
    while it <= n_max:
        x2 = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)))
        if abs(x2 - x1) < eps:
            return x2, it
        else:
            x0 = x1
            x1 = x2
        it += 1
    return False


def newton_raphson(f, x0, eps, n_max=1000):
    n = 1
    while n <= n_max:
        x1 = x0 - (f(x0) / derivative(f, x0))
        if abs(x1 - x0) < eps:
            return x1, n
        else:
            x0 = x1
        n += 1
    return False


start = process_time()
for i in range(0, 10000):
    root, iterations = bisect(function_1, -7, 7, 1e-5)
end = process_time()

print("For bisection method (x1 = {2}, x2={3}), elapsed time is: {0} s and number of iterations {1}".
      format(str((end - start) / 10000), str(iterations), str(-7), str(7)))
start = process_time()
for j in range(0, 10000):
    root, iterations = secant(function_1, -1, 7, 1e-5)
end = process_time()

print("For secant method,(x1 = {2}, x2={3}) elapsed time is: {0} s and number of iterations {1}".
      format(str((end - start) / 10000), str(iterations), str(-1), str(7)))

for k in range(0, 10000):
    root, iterations = newton_raphson(function_1, 4, 1e-5)
end = process_time()

print("For Newton_Raphson method (x1 ={2}, elapsed time is: {0} s and number of iterations {1}".
      format(str((end - start) / 10000), str(iterations), str(4)))

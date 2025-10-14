import numpy as np
from scipy.integrate import nquad
from scipy import stats
import random
from tqdm import tqdm
import time

def compute_etest(g, X, Y):
    X = np.asarray(X)
    Y = np.asarray(Y)

    # Создаем матрицы разностей
    diff_xy = X[:, np.newaxis] - Y[np.newaxis, :]  # форма: (len(X), len(Y))
    diff_xx = X[:, np.newaxis] - X[np.newaxis, :]  # форма: (len(X), len(X))
    diff_yy = Y[:, np.newaxis] - Y[np.newaxis, :]  # форма: (len(Y), len(Y))
    
    # Применяем функцию g ко всем элементам
    g_xy = g(diff_xy)
    g_xx = g(diff_xx)
    g_yy = g(diff_yy)

    # Суммируем все элементы (исключая диагональ для g_xx если нужно)
    phi_ab = np.sum(g_xy) 
    phi_a = np.sum(np.triu(g_xx, k=1)) 
    phi_b = np.sum(np.triu(g_yy, k=1)) 

    return (phi_ab - phi_a - phi_b) / X.shape[0]**2 


def compute_etest2(g, X, Y):
    n = len(X)

    phi_a = 0
    phi_b = 0
    phi_ab = 0

    for i in range(n):
        for j in range(n):
            if i < j:
                phi_a += g(X[i] - X[j])
                phi_b += g(Y[i] - Y[j])
            phi_ab += g(X[i] - Y[j])

    phi_a  /=  n**2
    phi_b  /=  n**2  #  /= m**2 in general case
    phi_ab /= n**2   #  /= n*m in general case

    phi = phi_ab - phi_a - phi_b

    return phi


def compute_integrals(f, g, d2_g):
    integrals = dict()

    print("computing J1...")
    integrals["J1"] = nquad(
        lambda x,y: g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J1 =", integrals["J1"])
    
    print("computing J2...")
    integrals["J2"] = nquad(
        lambda x,y: (g(x - y) ** 2) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J2 =", integrals["J2"])

    t0= time.time()
    print("computing J3...") # 22:51
    # integrals["J3"] = nquad(
    #     lambda x,y,z: g(x - y) * g(x - z) * f(x) * f(y) * f(z),
    #     ranges=[(-50, +50),
    #             (-50, +50),
    #             (-50, +50)], opts={"epsabs":1e-3, "epsrel":1e-3})[0]
    # integrals["J3"] = 0.763368 # for normal
    integrals["J3"] = 6.881056 # for cauchy
    print("J3 =", integrals["J3"])
    print("elapsed time for J3:", time.time()-t0)
    
    print("computing J1_star...")
    integrals["J1_star"] = 0.5 * nquad(
        lambda x,y: d2_g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J1_star =", integrals["J1_star"])

    print("computing J2_star...")
    integrals["J2_star"] = 0.5 * nquad(
        lambda x,y: (y**2 - 0.5 * (x-y)**2) * d2_g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])[0]
    print("J2_star =", integrals["J2_star"])


    return integrals



def compute_asymptotic_power(alpha, b, a):
    dist = stats.norm(loc=0.0, scale=1.0)
    F = dist.cdf
    z = dist.ppf(1 - alpha / 2)
    return 1 - F(z - b / a) + F( -z - b / a)


def compute_crit_val(n, M, alpha, template, g):
    
    dist1 = template(1,0,0)

    Z = dist1.rvs(2*n)

    test_vals = []
    indices = list(range(2*n))
    
    for _ in tqdm(range(M)):
        
        random.shuffle(indices)
        x_indices = indices[:n]
        y_indices = indices[n:]

        X = Z[x_indices]
        Y = Z[y_indices]
        etest_val = compute_etest(g, X, Y)
        test_vals.append( n * etest_val )

    return np.quantile(test_vals, 1 - alpha)


def compute_empirical_power(n, N, crit_val, d1, d2, g):
    """
    ``d1`` - first distribution \\
    ``d2`` - second distribution 
    """

    cnt_reject = 0

    for _ in tqdm(range(N)):
        X = d1.rvs(n)
        Y = d2.rvs(n)
        etest_val = compute_etest(g, X, Y)
        cnt_reject += int( n * etest_val >= crit_val)

    return cnt_reject / N 







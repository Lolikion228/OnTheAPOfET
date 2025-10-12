import numpy as np
from scipy.integrate import nquad
from scipy import stats

def compute_test(g, X, Y):
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


def compute_integrals(f, g, d2_g, h1, h2):
    integrals = dict()

    integrals["J1"] = nquad(
        lambda x,y: g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])
    
    integrals["J2"] = nquad(
        lambda x,y: (g(x - y) ** 2) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])
    
    integrals["J3"] = nquad(
        lambda x,y,z: g(x - y) * g(x - z) * f(x) * f(y) * f(z),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf),
                (-np.inf, +np.inf)])
    
    integrals["J1_star"] = 0.5 * (h1 ** 2) * nquad(
        lambda x,y: d2_g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])
    
    integrals["J2_star"] = 0.5 * (h2 ** 2) * nquad(
        lambda x,y: (y**2 - 0.5 * (x-y)**2) * d2_g(x - y) * f(x) * f(y),
        ranges=[(-np.inf, +np.inf),
                (-np.inf, +np.inf)])
    

    return integrals


def compute_asymptotic_power(alpha, b, a):
    pass


def compute_crit_val(n, M, alpha):
    pass


def compute_empirical_power(n, N, crit_val, d1, d2, g):
    """
    ``d1`` - first distribution \\
    ``d2`` - second distribution 
    """

    cnt_reject = 0

    for _ in range(N):
        X = d1.rvs(n)
        Y = d2.rvs(n)
        test_val = compute_test(g, X, Y)
        cnt_reject += int( n * test_val >= crit_val)

    return cnt_reject / N 







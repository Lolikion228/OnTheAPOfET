import numpy as np
from scipy.integrate import nquad


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


def compute_empirical_power(n, N, crit_val):
    pass







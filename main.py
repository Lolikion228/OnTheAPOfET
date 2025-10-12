from first import *



def g(x):
    pass

def d2_g(x):
    pass

# density for F(x)
def f(x):
    pass



np.random

def run_experiment(f, g, d2_g, 
                   h1=0, h2=2.1, alpha=0.05,
                   N=1000, M=17, sample_sizes=[100, 400, 900]):
    """
    ``h1`` - scale diff

    ``h2`` - shift diff

    ``alpha`` - level of significance 

    ``N`` - number of repeats computing empirical power

    ``M`` - number of permutations for computing critical value

    ``sample_sizes`` - sizes of samples to compute empirical power

    """
    integrals = compute_integrals(f, g, d2_g, h1, h2)

    b1 = ...
    b2 = ...
    
    a = ... 
    b = ...

    a_pow = compute_asymptotic_power(alpha, b, a)

    
    e_pows = []

    for n in sample_sizes:
        crit_val = compute_crit_val(n, M, alpha)
        e_pow = compute_empirical_power(n, N, crit_val)
        e_pows.append(e_pow)

    return e_pow, a_pow



def main():
    e_pow, a_pow = run_experiment(f=f, g=g, d2_g=d2_g)



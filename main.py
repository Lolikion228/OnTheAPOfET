from first import *
from scipy import stats


def g(x):
    pass

def d2_g(x):
    pass


#по имени и n генерировать выборки из F1, F2_n
# получать f(x) - плотность F1
def get_template(name):
    phi = None

    if name == "normal":
        def _phi(n, h1, h2):
            return stats.norm( ... ) 
        phi = _phi

    elif name == "cauchy":
        def _phi(n, h1, h2):
            return stats.cauchy( ... ) 
        phi = _phi
    
    return phi
    


# density for F(x)
def f(x):
    pass




def run_experiment(f, g, d2_g, 
                   h1=0, h2=2.1, alpha=0.05,
                   N=1000, M=17, sample_sizes=[100, 400, 900]):
    """
    params:
    ---
    ``h1`` - scale diff

    ``h2`` - shift diff

    ``alpha`` - level of significance 

    ``N`` - number of repeats computing empirical power

    ``M`` - number of permutations for computing critical value

    ``sample_sizes`` - sizes of samples to compute empirical power

    return:
    ---------
    ``empirical_powers`` and ``asymptotic_power``
    """

    integrals = compute_integrals(f, g, d2_g, h1, h2)

    b1 = np.sqrt( np.abs( integrals["J1_star"] ) )
    b2 = np.sqrt( np.abs( integrals["J2_star"] ) )
    
    a = (integrals["J2"] + integrals["J1"] ** 2 \
         - 2 * integrals["J3"]) ** 0.25
    b = np.sqrt( b1 ** 2 + b2 ** 2)

    emp_powers = []
    asp_power = compute_asymptotic_power(alpha, b, a)
    
    for n in sample_sizes:
        crit_val = compute_crit_val(n, M, alpha)
        e_pow = compute_empirical_power(n, N, crit_val)
        emp_powers.append(e_pow)

    return emp_powers, asp_power


def main():
    e_pow, a_pow = run_experiment(f=f, g=g, d2_g=d2_g)



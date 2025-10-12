from first import *
from scipy import stats
from tqdm import tqdm

def g(x):
    return np.log(1 + x**2)

def d2_g(x):
    return 2 * (1 - x**2) / ((1 + x**2)**2) 
 
templates = {
    "normal": lambda n,h1,h2: stats.norm(
                            loc = -h2 / (h1 + n**0.5),
                            scale= 1 / (1 + h1/(n**0.5)) ),
    "cauchy": None
}


def run_experiment(g, d2_g, template,
                   h1=0, h2=2.1, alpha=0.05,
                   N=1000, M=17, sample_sizes=[100, 400, 900]):
    """
    params:
    ---
    ``h1`` - scale diff (h1 >= 0)

    ``h2`` - shift diff

    ``alpha`` - level of significance 

    ``N`` - number of repeats computing empirical power

    ``M`` - number of permutations for computing critical value

    ``sample_sizes`` - sizes of samples to compute empirical power

    return:
    ---------
    ``empirical_powers`` and ``asymptotic_power``
    """

    if h1 < 0:
        raise Exception("h1 must be geq than 0")

    d1 = template(1, 0, 0)
    f = d1.pdf

    print("computing integrals...")
    integrals = compute_integrals(f, g, d2_g, h1, h2)
    print("done with integrals\n")

    b1 = np.sqrt( np.abs( integrals["J1_star"] ) )
    b2 = np.sqrt( np.abs( integrals["J2_star"] ) )
    
    a = (integrals["J2"] + integrals["J1"] ** 2 \
         - 2 * integrals["J3"]) ** 0.25
    b = np.sqrt( b1 ** 2 + b2 ** 2)

    emp_powers = []
    asp_power = compute_asymptotic_power(alpha, b, a)
    
    print("computing empirical_powers...")
    for n in sample_sizes:
        print(f"n={n}")
        print("computing critical_value...")
        crit_val = compute_crit_val(n, M, alpha, template, g)
        print("done computing critical_value")
        d2_n = template(n, h1, h2)
        print("computing e_pow for n*T_n...")
        e_pow = compute_empirical_power(n, N, crit_val, d1, d2_n, g)
        print("done computing e_pow for n*T_n\n")
        emp_powers.append(e_pow)

    return emp_powers, asp_power


def main():
    template = templates["normal"]
    e_pow, a_pow = run_experiment(g=g, d2_g=d2_g, template=template)
    print(e_pow)
    print(a_pow)

main()
import time

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Callable


def monte_carlo_with_delta(mc_integrand: Callable, delta_argument: Callable,
                           mc_bounds: Tuple[float, float],
                           n_bins: int, n_pts: int,
                           log=False):
    """Perform Monte Carlo integration when the integrand includes a delta function"""
    rand_pts_1 = np.random.uniform(0, 1, n_pts)
    rand_pts_2 = [np.random.uniform(1 - x1, 1) for x1 in rand_pts_1]
    rand_pts = np.array([rand_pts_1, rand_pts_2])
    rand_pts = rand_pts.T
    if log:
        bin_edges = np.geomspace(mc_bounds[0], mc_bounds[1], n_bins+1)
    else:
        bin_edges = np.linspace(mc_bounds[0], mc_bounds[1], n_bins+1)

    weights = np.array([mc_integrand(*pt) for pt in rand_pts])
    positions = np.array([delta_argument(*pt) for pt in rand_pts])

    # Sort the output arrays by the position in the output variable
    sort_indices = np.argsort(positions)
    weights = weights[sort_indices]
    positions = positions[sort_indices]

    bin_widths = bin_edges[1:] - bin_edges[:-1]

    # Bin the data, keeping track of errors as we go
    bin_totals = np.zeros(len(bin_edges)-1)
    uncertainty = np.zeros(len(bin_edges)-1)
    right_index = 0
    bin_index = 0
    while bin_index < len(bin_totals):
        if bin_index == 0:
            left_index = np.searchsorted(positions, bin_edges[0])
        else:
            left_index = right_index
        right_index = np.searchsorted(positions, bin_edges[bin_index+1])+1

        weights_in_bin = weights[left_index:right_index]
        bin_totals[bin_index] += sum(weights_in_bin)
        uncertainty[bin_index] += np.sqrt(sum(weights_in_bin ** 2))

        bin_index += 1

    return (bin_totals / bin_widths) / n_pts, (uncertainty / bin_widths) / n_pts, bin_edges


def heaviside_theta(arg):
    if arg > 0:
        return 1
    else:
        return 0


def integrand(x1, x2, z):
    x3 = 2 - x1 - x2
    min_val = min(x1, x2, x3)
    max_val = max(x1, x2, x3)

    ht1 = heaviside_theta(x1 + x2 - 1)
    ht2 = heaviside_theta(min_val / (2 - max_val) - z)
    if ht2 == 0:
        return 0

    main_portion = (x1 ** 2 + x2 ** 2) / ((1 - x1) * (1 - x2))
    return main_portion


def delta(x1, x2):
    x3 = 2 - x1 - x2
    max_val = max(x1, x2, x3)
    return 4 * (1 - max_val) / ((2 - max_val) ** 2)


def main():
    # z = 0.04
    n_bins = 40
    n_events = 500000
    for z in [0.04,
              #0.06, 0.08, 0.1
              ]:
        mc_pts, mc_unc, bin_edges = monte_carlo_with_delta(lambda x1, x2: integrand(x1, x2, z),
                                                           delta,
                                                           (0.005, 0.6), n_bins, n_events, log=True)
        x_pts = (bin_edges[:-1] + bin_edges[1:]) / 2

        plt.plot(x_pts, mc_pts * x_pts, label='z = %.2f' % z)
        plt.fill_between(x_pts, (mc_pts + mc_unc) * x_pts, (mc_pts - mc_unc) * x_pts, alpha=0.2)
    plt.xscale('log')
    plt.minorticks_on()
    plt.grid()
    plt.legend(fontsize=10)
    plt.ylabel(r'$C\,\rho\,\frac{d\sigma}{d\rho}$', fontsize=12)
    plt.xlabel(r'$\rho$', fontsize=12)
    plt.title((r'MC integration of $d\sigma / d\rho$' + '\nNo importance sampling'
               + '\n%i bins | %i samples' % (n_bins, n_events)), fontsize=12)
    plt.show()
    # plt.savefig('plots/mc_many_z_bins_%i_samples_%i.png' % (n_bins, n_events), bbox_inches='tight', padinches=0)
    # plt.savefig('plots/mc_many_z_bins_%i_samples_%i.pdf' % (n_bins, n_events), bbox_inches='tight', padinches=0)


if __name__ == '__main__':
    start = time.time()

    main()

    end = time.time()
    print('Executed in %s ms' % ((end - start) * 1000))

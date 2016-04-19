
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import getopt
import shelve


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "h", [])
        for opt, arg in opts:
            if opt == '-h':  # help
                print 'compare_results.py'
                sys.exit()
    except getopt.GetoptError:
        print 'compare_results.py'
        sys.exit(2)

    for a in (0, 1):
        for c in (0, 1):
            data_basic, gpr_basic, data_morew, dfp_morew = get_gpr_results(a, c)

            plot_slices(data_basic, gpr_basic, data_morew, dfp_morew)


def get_gpr_results(a, c):
    path = '../dw23_regression/sklearn_implementation/results/sk0.18/'

    filename = path + 'noise_d30w50_a%d_c%d.shelf' % (a, c)
    # filename = 'input/gpr_a%s_c%s.shelf' % (a, c)
    f = shelve.open(filename, flag='r')
    hist_basic = f['hist']
    dfp_basic = f['dfp']
    f.close()

    filename = path + 'd30w50_a%d_c%d.shelf' % (a, c)
    # filename = 'input/gpr_a%s_c%s.shelf' % (a, c)
    f = shelve.open(filename, flag='r')
    hist_morew = f['hist']
    hist_morew = hist_morew[(hist_morew['dw23_bin_center'] > -0.3) &
                            (hist_morew['dw23_bin_center'] < 0.3)]
    dfp_morew = f['dfp']
    f.close()

    return hist_basic, dfp_basic, hist_morew, dfp_morew


def plot_slices(hist1, dfp1, hist2, dfp2):
    plot_both = True
    # prepare dw23 slices in wness of the results for plotting
    data = hist1.groupby('wness_bin_center').agg(np.sum)['entries']*3

    gpr = dfp1.groupby('wness').agg(np.sum)['entries']

    gprs = dfp1.groupby('wness').agg(
        lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)

    gpr2 = dfp2.groupby('wness').agg(np.sum)['entries']

    gpr2s = dfp2.groupby('wness').agg(
        lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)

    # plot slices of the results
    fig, ax = plt.subplots()
    fig.set_figheight(38)
    fig.set_figwidth(15)
    ax.errorbar(data.index.values,
                data.values, np.sqrt(data).values,
                fmt='r.', markersize=10,
                elinewidth=2, label=u'Data')
    ax.plot(gpr, label='GPR w/ Noise')
    ax.fill(np.concatenate([gprs.index.values,
                            gprs.index.values[::-1]]),
            np.concatenate([gpr - 1.0 * gprs,
                            (gpr + 1.0 * gprs)[::-1]]),
            alpha=.5, fc='b', ec='None',
            label='Standard Deviation')
    if plot_both:
        ax.plot(gpr2, 'g', label='GPR')
        ax.fill(np.concatenate([gpr2s.index.values,
                                gpr2s.index.values[::-1]]),
                np.concatenate([gpr2 - 1.0 * gpr2s,
                                (gpr2 + 1.0 * gpr2s)[::-1]]),
                alpha=.5, fc='g', ec='None',
                label='Standard Deviation')
    ax.set_xlabel('Wness')
    ax.set_ylabel('Counts')
    ax.set_title('Fit Region Wness Distribution')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=16)

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])

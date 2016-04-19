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

            fcn = get_fcn_results(a, c, data_basic)

            plot_slices(data_basic, gpr_basic, data_morew, dfp_morew, fcn)


def get_gpr_results(a, c):
    path = '../dw23_regression/sklearn_implementation/results/sk0.18/'

    filename = path + 'd30w50_a%d_c%d.shelf' % (a, c)
    # filename = 'input/gpr_a%s_c%s.shelf' % (a, c)
    f = shelve.open(filename, flag='r')
    hist_basic = f['hist']
    dfp_basic = f['dfp']
    f.close()

    filename = path + 'morewing_d30w50_a%d_c%d.shelf' % (a, c)
    # filename = 'input/gpr_a%s_c%s.shelf' % (a, c)
    f = shelve.open(filename, flag='r')
    hist_morew = f['hist']
    hist_morew = hist_morew[(hist_morew['dw23_bin_center'] > -0.3) &
                            (hist_morew['dw23_bin_center'] < 0.3)]
    dfp_morew = f['dfp']
    f.close()

    return hist_basic, dfp_basic, hist_morew, dfp_morew


def read_fcn_fit_params(a, c, method):
    filename = 'input/2d_fit_params.csv'
    df = pd.read_csv(filename)
    return df[(df['arm'] == a) & (df['charge'] == c) & (df['method'] == method)]


def eval_2d_gaus_offset(par, d, w):
    gaus_mean = par[0] + par[1]*w
    mean_offset = par[13] + par[14]*w
    gaus_sigma1 = par[2] + par[3]*w  # should be wider
    gaus_sigma2 = par[4] + par[5]*w  # should be narrower
    gaus_factor = par[6] + par[7]*w
    pi = 3.14159265358979323846

    # polynomial describes the normalized 1D wness distribution.
    # used as a weight factor for the otherwise normalized dw23
    # distribution and different regions of wness
    polynomial_factor = par[8] + par[9]*w + par[10]*w**2 + \
        par[11]*w**3 + par[12]*w**4
    normalize_factor = \
        1 / (gaus_sigma1*np.sqrt(2*pi) + gaus_factor*gaus_sigma2*np.sqrt(2*pi))

    val = polynomial_factor * normalize_factor * \
        (
         np.exp(-0.5*((d-gaus_mean)/gaus_sigma1)**2) +
         gaus_factor*np.exp(-0.5*((d-(gaus_mean+mean_offset))/gaus_sigma2)**2)
         )

    return val


def get_fcn_results(a, c, data):
    # a = 1 - a
    par_df = read_fcn_fit_params(a, c, 'offset')

    pars = {}
    for n, row in par_df.iterrows():
        pars[row['param_index']] = row['param_value']

    d = []
    w = []
    e = []
    for n, row in data.iterrows():
        d.append(row['dw23_bin_center'])
        w.append(row['wness_bin_center'])
        e.append(eval_2d_gaus_offset(pars, d[-1], w[-1]))

    df = pd.DataFrame({'dw23': d, 'wness': w, 'entries': e})

    # scale the fcn data by the integral of the data (it was generated
    # from a normalized histogram of the data)
    # NOTE: these values taken from the c code where the fits were done
    if a == 0:
        if c == 0:
            integral = 63.942402
        else:
            integral = 69.702400
    else:
        if c == 0:
            integral = 89.776802
        else:
            integral = 94.155197
    # factor of 2 to account for 60 dw23 binning vs 30 there
    df['entries'] = df['entries'] * integral

    return df


def plot_slices(data, gpr, data_mw, gpr_mw, fcn, plot_both=False):
    # prepare dw23 slices in wness of the results for plotting
    data2 = data[(data['wness_bin_center'] > .1) &
                 (data['wness_bin_center'] < .3)]\
        .groupby('dw23_bin_center').agg(np.sum)['entries']
    data4 = data[(data['wness_bin_center'] > .3) &
                 (data['wness_bin_center'] < .5)]\
        .groupby('dw23_bin_center').agg(np.sum)['entries']
    data6 = data[(data['wness_bin_center'] > .5) &
                 (data['wness_bin_center'] < .7)]\
        .groupby('dw23_bin_center').agg(np.sum)['entries']
    data8 = data[(data['wness_bin_center'] > .7) &
                 (data['wness_bin_center'] < .9)]\
        .groupby('dw23_bin_center').agg(np.sum)['entries']

    gpr2 = gpr[(gpr['wness'] > .1) & (gpr['wness'] < .3)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr4 = gpr[(gpr['wness'] > .3) & (gpr['wness'] < .5)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr6 = gpr[(gpr['wness'] > .5) & (gpr['wness'] < .7)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr8 = gpr[(gpr['wness'] > .7) & (gpr['wness'] < .9)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr9 = gpr[(gpr['wness'] > .92) & (gpr['wness'] < 1.0)]\
        .groupby('dw23').agg(np.sum)['entries']

    gpr2s = gpr[(gpr['wness'] > .1) & (gpr['wness'] < .3)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr4s = gpr[(gpr['wness'] > .3) & (gpr['wness'] < .5)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr6s = gpr[(gpr['wness'] > .5) & (gpr['wness'] < .7)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr8s = gpr[(gpr['wness'] > .7) & (gpr['wness'] < .9)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr9s = gpr[(gpr['wness'] > .92) & (gpr['wness'] < 1.0)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)

    gpr_mw2 = gpr_mw[(gpr_mw['wness'] > .1) & (gpr_mw['wness'] < .3)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr_mw4 = gpr_mw[(gpr_mw['wness'] > .3) & (gpr_mw['wness'] < .5)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr_mw6 = gpr_mw[(gpr_mw['wness'] > .5) & (gpr_mw['wness'] < .7)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr_mw8 = gpr_mw[(gpr_mw['wness'] > .7) & (gpr_mw['wness'] < .9)]\
        .groupby('dw23').agg(np.sum)['entries']
    gpr_mw9 = gpr_mw[(gpr_mw['wness'] > .92) & (gpr_mw['wness'] < 1.0)]\
        .groupby('dw23').agg(np.sum)['entries']

    gpr_mw2s = gpr_mw[(gpr_mw['wness'] > .1) & (gpr_mw['wness'] < .3)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr_mw4s = gpr_mw[(gpr_mw['wness'] > .3) & (gpr_mw['wness'] < .5)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr_mw6s = gpr_mw[(gpr_mw['wness'] > .5) & (gpr_mw['wness'] < .7)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr_mw8s = gpr_mw[(gpr_mw['wness'] > .7) & (gpr_mw['wness'] < .9)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
    gpr_mw9s = gpr_mw[(gpr_mw['wness'] > .92) & (gpr_mw['wness'] < 1.0)]\
        .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)

    fcn2 = fcn[(fcn['wness'] > .1) & (fcn['wness'] < .3)]\
        .groupby('dw23').agg(np.sum)['entries']
    fcn4 = fcn[(fcn['wness'] > .3) & (fcn['wness'] < .5)]\
        .groupby('dw23').agg(np.sum)['entries']
    fcn6 = fcn[(fcn['wness'] > .5) & (fcn['wness'] < .7)]\
        .groupby('dw23').agg(np.sum)['entries']
    fcn8 = fcn[(fcn['wness'] > .7) & (fcn['wness'] < .9)]\
        .groupby('dw23').agg(np.sum)['entries']

    # plot slices of the results
    fig, ax = plt.subplots(2, 3)
    fig.set_figheight(38)
    fig.set_figwidth(15)
    ax[0, 0].errorbar(data2.index.values,
                      data2.values, np.sqrt(data2).values,
                      fmt='r.', markersize=10,
                      elinewidth=2, label=u'Observations')
    ax[0, 0].plot(gpr2)
    ax[0, 0].fill(np.concatenate([gpr2s.index.values,
                                  gpr2s.index.values[::-1]]),
                  np.concatenate([gpr2 - 1.9600 * gpr2s,
                                  (gpr2 + 1.9600 * gpr2s)[::-1]]),
                  alpha=.5, fc='b', ec='None',
                  label='95% confidence interval')
    if plot_both:
        ax[0, 0].plot(gpr_mw2, 'g')
        ax[0, 0].fill(np.concatenate([gpr_mw2s.index.values,
                                      gpr_mw2s.index.values[::-1]]),
                      np.concatenate([gpr_mw2 - 1.9600 * gpr_mw2s,
                                      (gpr_mw2 + 1.9600 * gpr_mw2s)[::-1]]),
                      alpha=.5, fc='g', ec='None',
                      label='95% confidence interval')
    ax[0, 0].plot(fcn2, 'black')

    ax[0, 1].errorbar(data4.index.values,
                      data4.values, np.sqrt(data4).values,
                      fmt='r.', markersize=10,
                      elinewidth=2, label=u'Observations')
    ax[0, 1].plot(gpr4)
    ax[0, 1].fill(np.concatenate([gpr4s.index.values,
                                  gpr4s.index.values[::-1]]),
                  np.concatenate([gpr4 - 1.9600 * gpr4s,
                                  (gpr4 + 1.9600 * gpr4s)[::-1]]),
                  alpha=.5, fc='b', ec='None',
                  label='95% confidence interval')
    if plot_both:
        ax[0, 1].plot(gpr_mw4, 'g')
        ax[0, 1].fill(np.concatenate([gpr_mw4s.index.values,
                                      gpr_mw4s.index.values[::-1]]),
                      np.concatenate([gpr_mw4 - 1.9600 * gpr_mw4s,
                                      (gpr_mw4 + 1.9600 * gpr_mw4s)[::-1]]),
                      alpha=.5, fc='g', ec='None',
                      label='95% confidence interval')
    ax[0, 1].plot(fcn4, 'black')

    ax[0, 2].errorbar(data6.index.values,
                      data6.values, np.sqrt(data6).values,
                      fmt='r.', markersize=10,
                      elinewidth=2, label=u'Observations')
    ax[0, 2].plot(gpr6)
    ax[0, 2].fill(np.concatenate([gpr6s.index.values,
                                  gpr6s.index.values[::-1]]),
                  np.concatenate([gpr6 - 1.9600 * gpr6s,
                                  (gpr6 + 1.9600 * gpr6s)[::-1]]),
                  alpha=.5, fc='b', ec='None',
                  label='95% confidence interval')
    if plot_both:
        ax[0, 2].plot(gpr_mw6, 'g')
        ax[0, 2].fill(np.concatenate([gpr_mw6s.index.values,
                                      gpr_mw6s.index.values[::-1]]),
                      np.concatenate([gpr_mw6 - 1.9600 * gpr_mw6s,
                                      (gpr_mw6 + 1.9600 * gpr_mw6s)[::-1]]),
                      alpha=.5, fc='g', ec='None',
                      label='95% confidence interval')
    ax[0, 2].plot(fcn6, 'black')

    ax[1, 0].errorbar(data8.index.values,
                      data8.values, np.sqrt(data8).values,
                      fmt='r.', markersize=10,
                      elinewidth=2, label=u'Observations')
    ax[1, 0].plot(gpr8)
    ax[1, 0].fill(np.concatenate([gpr8s.index.values,
                                  gpr8s.index.values[::-1]]),
                  np.concatenate([gpr8 - 1.9600 * gpr8s,
                                  (gpr8 + 1.9600 * gpr8s)[::-1]]),
                  alpha=.5, fc='b', ec='None',
                  label='95% confidence interval')
    if plot_both:
        ax[1, 0].plot(gpr_mw8, 'g')
        ax[1, 0].fill(np.concatenate([gpr_mw8s.index.values,
                                      gpr_mw8s.index.values[::-1]]),
                      np.concatenate([gpr_mw8 - 1.9600 * gpr_mw8s,
                                      (gpr_mw8 + 1.9600 * gpr_mw8s)[::-1]]),
                      alpha=.5, fc='g', ec='None',
                      label='95% confidence interval')
    ax[1, 0].plot(fcn8, 'black')

    ax[1, 1].plot(gpr9)
    ax[1, 1].fill(np.concatenate([gpr9s.index.values,
                                  gpr9s.index.values[::-1]]),
                  np.concatenate([gpr9 - 1.9600 * gpr9s,
                                  (gpr9 + 1.9600 * gpr9s)[::-1]]),
                  alpha=.5, fc='b', ec='None',
                  label='95% confidence interval')
    if plot_both:
        ax[1, 1].plot(gpr_mw9, 'g')
        ax[1, 1].fill(np.concatenate([gpr_mw9s.index.values,
                                      gpr_mw9s.index.values[::-1]]),
                      np.concatenate([gpr_mw9 - 1.9600 * gpr_mw9s,
                                      (gpr_mw9 + 1.9600 * gpr_mw9s)[::-1]]),
                      alpha=.5, fc='g', ec='None',
                      label='95% confidence interval')

    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])

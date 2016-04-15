import shelve
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

np.random.seed(1)


class data_instance:
    '''
    this class object contains a subset of data and functions
    to process this data
    '''

    def __init__(self, a, c,
                 d_nbins=60, w_nbins=50,
                 narrow_fit=False,
                 flatten_wness=False,
                 bigwidebins=False,
                 morewings=False):
        self.arm = a
        self.charge = c
        self.d_nbins = d_nbins
        self.w_nbins = w_nbins
        self.narrow_fit = narrow_fit
        self.flatten_wness = flatten_wness
        self.bigwidebins = bigwidebins
        self.morewings = morewings

        tags = ''
        if narrow_fit:
            tags = tags + 'narrow_'
        if flatten_wness:
            tags = tags + 'flatw_'
        if bigwidebins:
            tags = tags + 'bwb_'
        if morewings:
            tags = tags + 'morewing_'

        self.filename = 'results/sk0.18/%sd%dw%d_a%s_c%s.shelf' \
            % (tags, d_nbins, w_nbins, a, c)

    def set_data(self, df):
        self.df = df

    def build_dw23wness_hist(self):
        # create a binned histogram of entries with respect to dw23 and wness
        dw23bins = np.linspace(-.3, .3, self.d_nbins+1)
        wnessbins = np.linspace(0, 1, self.w_nbins+1)
        self.hist = self.df.groupby([pd.cut(self.df['dw23'], dw23bins),
                                     pd.cut(self.df['Wness'], wnessbins)]
                                    ).count()
        self.hist['entries'] = self.hist['Wness']
        self.hist = self.hist[['entries', 'pz']]
        self.hist.reset_index(inplace=True)
        self.hist.drop('pz', axis=1, inplace=True)
        self.hist.fillna(0, inplace=True)

        def convert_range_str(s):
            '''helper function:
            converts the text discribing bin range to a bin center value'''
            s = s[1:-1]
            low, up = s.split(', ')

            mid = (float(low) + float(up)) / 2

            return mid

        self.hist['dw23_bin_center'] = self.hist['dw23'].apply(
            convert_range_str)
        self.hist['wness_bin_center'] = self.hist['Wness'].apply(
            convert_range_str)

        # trim the extremes of the wness range
        self.hist = self.hist[(self.hist['wness_bin_center'] > .1) &
                              (self.hist['wness_bin_center'] < .9)]
        if self.narrow_fit:
            self.hist = self.hist[(self.hist['dw23_bin_center'] > -.1) &
                                  (self.hist['dw23_bin_center'] < .1)]

        print hist.

    def do_flatten_wness(self):
        print 'success!'
        self.slice_sums = self.hist.groupby('wness_bin_center').agg(np.sum)

        y = []
        for n, row in self.hist.iterrows():
            scaled_value = \
                row['entries'] / \
                self.slice_sums.loc[row['wness_bin_center']]['entries']

            y.append(scaled_value)

        return y

    def fit_gpr_model(self):
        '''
        create and fit the gaussian process model.
        the results in the model object are stored in a class variable
        '''
        # prepare input for fitting function
        X = self.hist[['wness_bin_center', 'dw23_bin_center']].values
        y = self.hist['entries'].values

        # uncertainty of counts from a poisson distribution
        dy = np.sqrt(y)
        # "nugget" is used to inform fitting algorithm of input uncertainty
        nugget = (dy / y) ** 2
        inds = np.where(np.isnan(nugget))
        nugget[inds] = 1.0

        if self.flatten_wness:
            y = self.do_flatten_wness()

        # define kernel
        self.kernel = RBF([.1, .1])

        # Instanciate a Gaussian Process model
        self.gp = GaussianProcessRegressor(kernel=self.kernel,
                                           alpha=nugget,
                                           normalize_y=True,
                                           n_restarts_optimizer=50)

        # Fit to data using Maximum Likelihood Estimation of the parameters
        self.gp.fit(X, y)

    def do_unflatten_wness(self, x, y_pred, sigma):
        print type(x), x.shape
        print type(y_pred), y_pred.shape
        y_unflat = []
        for i, coords in enumerate(x):
            print coords[0]
            if coords[0] in self.slice_sums.index:
                y_unflat.append(y_pred[i]*self.slice_sums[coords[0]]['entries'])
            else:
                done = False
                next_w = coords[0]
                while not done:
                    next_w += .001
                    if next_w in self.slice_sums.index:
                        done = True
                    # if next_w > max ---

        return y_pred, sigma

    def make_predictions(self):
        '''
        use the fitted model to predict results
        '''
        # change
        # d_binsize = 0.6 / self.d_nbins
        d_binsize = 0.6 / 90
        w_binsize = 1.0 / self.w_nbins
        # generate grid of points for which to make predicitions
        x = []
        for w in np.arange(.1+w_binsize/2, 1.0+w_binsize/2, w_binsize):
            for d in np.arange(-.3+d_binsize/2, .3+d_binsize/2, d_binsize):
                x.append([w, d])
        x = np.array(x)

        # make predictions
        y_pred, sigma = self.gp.predict(x, return_std=True)

        if self.flatten_wness:
            y_pred, sigma = self.do_unflatten_wness(x, y_pred, sigma)

        # create dataframe of prediction results
        self.dfp = pd.DataFrame({'wness': x[:, 0],
                                 'dw23': x[:, 1],
                                 'entries': y_pred,
                                 'sigma': sigma})

    def get_chi2ndf(self, dw23range='full'):
        '''
        calculate and return the chi2 / ndf of the fit results
        relative to the data points
        '''
        chi2 = 0
        count = 0
        rev_chi2 = 0
        rev_count = 0

        # cast these values as strings to avoid faulty float comparison
        # they will be returned to float at the end
        self.dfp[['wness', 'dw23']] = self.dfp[['wness', 'dw23']].applymap(str)

        for n, row in self.hist.iterrows():
            dd = row['dw23_bin_center']  # dw23 from data

            if dw23range == 'full' or\
                    (dw23range == 'narrow' and dd > -.1 and dd < .1) or\
                    (dw23range == 'wings' and (dd < -.1 or dd > .1)):
                dd = str(dd)  # dw23 from data
                wd = str(row['wness_bin_center'])  # wness from data
                ed = row['entries']  # entries from data
                ud = np.sqrt(ed)  # uncertainty from data
                if ud == 0:
                    ud = 1

                # entries from fit
                ef = self.dfp[(self.dfp['dw23'] == dd) &
                              (self.dfp['wness'] == wd)][
                                  'entries'].values[0]
                uf = self.dfp[(self.dfp['dw23'] == dd) &
                              (self.dfp['wness'] == wd)][
                                    'sigma'].values[0]

                try:

                    chi2 += (ed - ef)**2 / ud**2
                    count += 1
                except:
                    print 'problem with dw23=%.4f wness=%.4f count%d'\
                        % (dd, wd, count)

                try:
                    rev_chi2 += (ef - ed)**2 / uf**2
                    rev_count += 1
                except:
                    print 'problem with dw23=%.4f wness=%.4f count%d'\
                        % (dd, wd, rev_count)

        # change back the values from strings to floats
        self.dfp[['wness', 'dw23']] = \
            self.dfp[['wness', 'dw23']].applymap(float)

        return chi2/(count-1), rev_chi2/(rev_count-1)

    def save_results(self):
        f = shelve.open(self.filename, flag='n')
        f['hist'] = self.hist
        f['dfp'] = self.dfp
        f['gp'] = self.gp
        f['kernel'] = self.kernel
        f.close()

    def export_results(self):
        filename = 'results/gpr_predictions_a%sc%s.csv' % (
            self.arm, self.charge)
        self.dfp.to_csv(filename, index=False)

    def read_results(self):
        f = shelve.open(self.filename, flag='r')
        self.hist = f['hist']
        self.dfp = f['dfp']
        f.close()

    def plot_2d_hist(self):
        # create lists for plotting
        wness = self.hist['wness_bin_center'].values
        dw23 = self.hist['dw23_bin_center'].values
        entries = self.hist['entries'].values

        # plot histogram
        plt.figure(num=1, figsize=(9, 7))
        ax = plt.subplot(projection='3d')
        ax.scatter(wness, dw23, entries)
        ax.view_init(elev=30, azim=-160)
        ax.set_xlabel('wness')
        ax.set_ylabel('dw23')
        ax.set_zlabel('counts')
        plt.show()

    def plot_raw_dw23_slices(self):
        # select subsets of the dw23 range for plotting slices
        df2 = self.hist[(self.hist['wness_bin_center'] > .1) &
                        (self.hist['wness_bin_center'] < .3)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df4 = self.hist[(self.hist['wness_bin_center'] > .3) &
                        (self.hist['wness_bin_center'] < .5)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df6 = self.hist[(self.hist['wness_bin_center'] > .5) &
                        (self.hist['wness_bin_center'] < .7)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df8 = self.hist[(self.hist['wness_bin_center'] > .7) &
                        (self.hist['wness_bin_center'] < .9)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']

        fig, ax = plt.subplots(2, 2)
        fig.set_figheight(10)
        fig.set_figwidth(15)
        ax[0, 0].plot(df2)
        ax[0, 1].plot(df4)
        ax[1, 0].plot(df6)
        ax[1, 1].plot(df8)
        plt.show()

    def generate_pct_err_curve(self, df2p, df2ps, df4p, df4ps,
                               df6p, df6ps, df8p, df8ps):
        pct_err2 = np.absolute(df2ps.values / df2p.values)
        pct_err4 = np.absolute(df4ps.values / df4p.values)
        pct_err6 = np.absolute(df6ps.values / df6p.values)
        pct_err8 = np.absolute(df8ps.values / df8p.values)

        avg_pct_err = (pct_err2 + pct_err4 + pct_err6 + pct_err8) / 4.0 * 100.0

        return pd.Series(avg_pct_err, index=df2p.index)

    def plot_result_dw23_slices(self):
        if self.narrow_fit:
            self.dfp = self.dfp[(self.dfp['dw23'] > -.1) &
                                (self.dfp['dw23'] < .1)]

        # prepare dw23 slices in wness of the results for plotting
        df2 = self.hist[(self.hist['wness_bin_center'] > .1) &
                        (self.hist['wness_bin_center'] < .3)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df4 = self.hist[(self.hist['wness_bin_center'] > .3) &
                        (self.hist['wness_bin_center'] < .5)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df6 = self.hist[(self.hist['wness_bin_center'] > .5) &
                        (self.hist['wness_bin_center'] < .7)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']
        df8 = self.hist[(self.hist['wness_bin_center'] > .7) &
                        (self.hist['wness_bin_center'] < .9)]\
            .groupby('dw23_bin_center').agg(np.sum)['entries']

        df2p = self.dfp[(self.dfp['wness'] > .1) & (self.dfp['wness'] < .3)]\
            .groupby('dw23').agg(np.sum)['entries']
        df4p = self.dfp[(self.dfp['wness'] > .3) & (self.dfp['wness'] < .5)]\
            .groupby('dw23').agg(np.sum)['entries']
        df6p = self.dfp[(self.dfp['wness'] > .5) & (self.dfp['wness'] < .7)]\
            .groupby('dw23').agg(np.sum)['entries']
        df8p = self.dfp[(self.dfp['wness'] > .7) & (self.dfp['wness'] < .9)]\
            .groupby('dw23').agg(np.sum)['entries']
        df9p = self.dfp[(self.dfp['wness'] > .9) & (self.dfp['wness'] < 1)]\
            .groupby('dw23').agg(np.sum)['entries']

        df2ps = self.dfp[(self.dfp['wness'] > .1) & (self.dfp['wness'] < .3)]\
            .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
        df4ps = self.dfp[(self.dfp['wness'] > .3) & (self.dfp['wness'] < .5)]\
            .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
        df6ps = self.dfp[(self.dfp['wness'] > .5) & (self.dfp['wness'] < .7)]\
            .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
        df8ps = self.dfp[(self.dfp['wness'] > .7) & (self.dfp['wness'] < .9)]\
            .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)
        df9ps = self.dfp[(self.dfp['wness'] > .9) & (self.dfp['wness'] < 1)]\
            .groupby('dw23').agg(lambda x: np.sum(x**2))['sigma'].apply(np.sqrt)

        # plot slices of the results
        fig, ax = plt.subplots(2, 3)
        fig.set_figheight(38)
        fig.set_figwidth(15)
        ax[0, 0].errorbar(df2.index.values,
                          df2.values, np.sqrt(df2).values,
                          fmt='r.', markersize=10,
                          elinewidth=2, label=u'Observations')
        ax[0, 0].plot(df2p)
        ax[0, 0].fill(np.concatenate([df2ps.index.values,
                                      df2ps.index.values[::-1]]),
                      np.concatenate([df2p - 1.9600 * df2ps,
                                      (df2p + 1.9600 * df2ps)[::-1]]),
                      alpha=.5, fc='b', ec='None',
                      label='95% confidence interval')
        ax[0, 1].errorbar(df4.index.values,
                          df4.values, np.sqrt(df4).values,
                          fmt='r.', markersize=10,
                          elinewidth=2, label=u'Observations')
        ax[0, 1].plot(df4p)
        ax[0, 1].fill(np.concatenate([df4ps.index.values,
                                      df4ps.index.values[::-1]]),
                      np.concatenate([df4p - 1.9600 * df4ps,
                                      (df4p + 1.9600 * df4ps)[::-1]]),
                      alpha=.5, fc='b', ec='None',
                      label='95% confidence interval')
        ax[0, 2].errorbar(df6.index.values,
                          df6.values, np.sqrt(df6).values,
                          fmt='r.', markersize=10,
                          elinewidth=2, label=u'Observations')
        ax[0, 2].plot(df6p)
        ax[0, 2].fill(np.concatenate([df6ps.index.values,
                                      df6ps.index.values[::-1]]),
                      np.concatenate([df6p - 1.9600 * df6ps,
                                      (df6p + 1.9600 * df6ps)[::-1]]),
                      alpha=.5, fc='b', ec='None',
                      label='95% confidence interval')
        ax[1, 0].errorbar(df8.index.values,
                          df8.values, np.sqrt(df8).values,
                          fmt='r.', markersize=10,
                          elinewidth=2, label=u'Observations')
        ax[1, 0].plot(df8p)
        ax[1, 0].fill(np.concatenate([df8ps.index.values,
                                      df8ps.index.values[::-1]]),
                      np.concatenate([df8p - 1.9600 * df8ps,
                                      (df8p + 1.9600 * df8ps)[::-1]]),
                      alpha=.5, fc='b', ec='None',
                      label='95% confidence interval')
        ax[1, 1].plot(df9p)
        ax[1, 1].fill(np.concatenate([df9ps.index.values,
                                      df9ps.index.values[::-1]]),
                      np.concatenate([df9p - 1.9600 * df9ps,
                                      (df9p + 1.9600 * df9ps)[::-1]]),
                      alpha=.5, fc='b', ec='None',
                      label='95% confidence interval')

        # generate and plot percent error curve
        dferr = self.generate_pct_err_curve(df2p, df2ps, df4p, df4ps,
                                            df6p, df6ps, df8p, df8ps)
        ax[1, 2].plot(dferr)
        ax[1, 2].set_yscale('log')

        plt.show()

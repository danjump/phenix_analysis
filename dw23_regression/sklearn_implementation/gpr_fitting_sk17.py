import sys
import getopt
import pandas as pd

from data_instance import data_instance


def split_ac(df):
    '''
    input is dataframe of raw input data
    output is a 2x2 list of dataframes separated by arm/charge
    '''
    # adjust arm and charge values to 0/1
    df['charge'] = (df['charge'] > 0).apply(int)
    df['arm'] = (df['pz'] > 0).apply(int)

    # create a list of dfs separated for each arm/charge
    grouped = df.groupby(['arm', 'charge'])
    dfs = [[None for x in range(2)] for x in range(2)]
    for g in grouped:
        dfs[g[0][0]][g[0][1]] = g[1]

    return dfs


def main(argv):
    # --------------------------------------------------------------------------
    #                    PARSE INPUT ARGUMENTS
    # --------------------------------------------------------------------------

    # intitialize default arguments
    which_arms = ['0', '1']
    which_charge = ['0', '1']
    skipfit = False
    simulate_error = False
    plot = False
    narrow_fit = False
    d_nbins = 60

    # get/parse command line arguments and set variables accordingly
    try:
        opts, args = getopt.getopt(argv, "fsb:pnha:c:",
                                   ['verbose', 'plot', 'dnbins=',
                                    'narrow', 'skipfit', 'simerr'])
        for opt, arg in opts:
            if opt == '-h':  # help
                print 'gpr_fitting.py -a <armoption> -c <chargeoption>'
                sys.exit()
            elif opt == '-a':  # arm option
                # the expected input argument is 0, 1 or 2 for south
                # north or both respectively
                # which_arms will end up being a list of characters
                # containing one or both of '0' or '1'
                if int(arg) < 2:
                    which_arms = [arg]
                elif int(arg) == 2:
                    which_arms = ['0', '1']
                else:
                    print 'ERROR: bad arm option!'
                    sys.exit()
            elif opt == '-c':  # charge option
                # the expected input argument is 0, 1 or 2 for negative
                # positive or both respectively
                # which_charge will end up being a list of characters
                # containing one or both of '0' or '1'
                if int(arg) < 2:
                    which_charge = [arg]
                elif int(arg) > 2:
                    which_charge = ['0', '1']
                else:
                    print 'ERROR: bad charge option!'
                    sys.exit()
            elif opt in ('-f', '--skipfit'):
                skipfit = True  # read previous results of fit
            elif opt in ('-s', '--simerr'):
                simulate_error = True  # read previous results of fit
            elif opt in ('-p', '--plot'):
                plot = True  # do plotting
            elif opt in ('-n', '--narrow'):
                narrow_fit = True  # use -.1 < dw23 < .1
            elif opt in ('-b', '--dnbins'):
                d_nbins = int(arg)  # number of dw23 bins
    except getopt.GetoptError:
        print 'test.py -a <armoption> -c <chargeoption>'
        sys.exit(2)

    # ---------------------------------------------------------------------
    #                         LOAD DATA
    # ---------------------------------------------------------------------
    if not skipfit:
        # read data
        df = pd.read_csv('input/data_wness.csv', index_col='index')

        # segment to a list of dfs separated by a/c
        dfs = split_ac(df)

    # master loop for arms/charges
    for a in which_arms:
        for c in which_charge:
            # initialize data objects
            data_obj = data_instance(a, c,
                                     d_nbins=d_nbins,
                                     narrow_fit=narrow_fit,
                                     simulate_error=simulate_error)

            if not skipfit:
                data_obj.set_data(dfs[int(a)][int(c)])

            # -----------------------------------------------------------
            #                  DO FITTING
            # -----------------------------------------------------------
            if skipfit:
                data_obj.read_results()
            else:
                data_obj.build_dw23wness_hist()
                data_obj.fit_gpr_model()
                data_obj.make_predictions()
                data_obj.save_results()

            chi2, rev_chi2 = data_obj.get_chi2ndf()
            nar_chi2, nar_rev_chi2 = data_obj.get_chi2ndf('narrow')
            wing_chi2, wing_rev_chi2 = data_obj.get_chi2ndf('wings')
            print 'ac=%s chi2ndf=%.4f narchi2ndf=%.4f wingchi2ndf=%.4f' % \
                (a+c, chi2, nar_chi2, wing_chi2)
            print 'reverse ac=%s chi2ndf=%.4f narchi2ndf=%.4f wingchi2ndf=%.4f'\
                % (a+c, rev_chi2, nar_rev_chi2, wing_rev_chi2)

            # -----------------------------
            #          PLOTTING
            # -----------------------------
            if plot:
                data_obj.plot_result_dw23_slices()

            data_obj.export_results()


    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv[1:])

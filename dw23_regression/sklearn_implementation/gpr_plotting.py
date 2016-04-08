import sys
import getopt
import pandas as pd

from data_instance import data_instance


def main(argv):
    # --------------------------------------------------------------------------
    #                    PARSE INPUT ARGUMENTS
    # --------------------------------------------------------------------------

    # intitialize default arguments
    which_arms = ['0', '1']
    which_charge = ['0', '1']
    simulate_error = False
    narrow_fit = False
    d_nbins = 60

    # get/parse command line arguments and set variables accordingly
    try:
        opts, args = getopt.getopt(argv, "sb:nha:c:",
                                   ['verbose', 'dnbins=',
                                    'narrow', 'simerr'])
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
            elif opt in ('-s', '--simerr'):
                simulate_error = True  # read previous results of fit
            elif opt in ('-n', '--narrow'):
                narrow_fit = True  # use -.1 < dw23 < .1
            elif opt in ('-b', '--dnbins'):
                d_nbins = int(arg)  # number of dw23 bins
    except getopt.GetoptError:
        print 'test.py -a <armoption> -c <chargeoption>'
        sys.exit(2)
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow_simerr-10.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow_simerr.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_simerr.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d60w50.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c1_d60w50.shelf/home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow_simerr-10.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_narrow_simerr.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d30w50_simerr.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c0_d60w50.shelf
        /home/danielj/work/phenix_analysis/run13_analysis/dw23_regression/sklearn_implementation/results/a0_c1_d60w50.shelf

    for a in which_arms:
        for c in which_charge:
            # initialize data objects
            data_obj = data_instance(a, c,
                                     d_nbins=d_nbins,
                                     narrow_fit=narrow_fit,
                                     simulate_error=simulate_error)

            data_obj.read_results()

            chi2, rev_chi2 = data_obj.get_chi2ndf()
            nar_chi2, nar_rev_chi2 = data_obj.get_chi2ndf('narrow')
            wing_chi2, wing_rev_chi2 = data_obj.get_chi2ndf('wings')
            print 'ac=%s chi2ndf=%.4f narchi2ndf=%.4f wingchi2ndf=%.4f' % \
                (a+c, chi2, nar_chi2, wing_chi2)
            print 'reverse ac=%s chi2ndf=%.4f narchi2ndf=%.4f wingchi2ndf=%.4f'\
                % (a+c, rev_chi2, nar_rev_chi2, wing_rev_chi2)

            data_obj.plot_result_dw23_slices()

    sys.exit(0)

if __name__ == "__main__":
    main(sys.argv[1:])

import numpy as np
import pandas as pd

from gpr_instance import gpr_instance

class gpr_analysis_obj:
    """a master class for the grp analysis

    contains functions and variables to facilitate
    the flow of executing the analysisi
    """

    def __init__(self, args_dict):
        """initialize class object with analysis options"""
        self.do_method_raw =    args_dict['do_method_raw'   ]
        self.do_method_scaled = args_dict['do_method_scaled']
        self.do_method_log =    args_dict['do_method_log'   ]
        self.do_range_narrow =  args_dict['do_range_narrow' ]
        self.do_range_full =    args_dict['do_range_full'   ]
        self.do_threshold =     args_dict['do_threshold'    ]
        self.input_file =       args_dict['input_file'      ]

        # generate array of thresholds that we want to do (True ones only)
        self.true_thrs = self.do_threshold[self.do_threshold[:, 1]==True][:, 0]

        # zero-initialize arrays that will hold gpr_instance objects
        # 3 indecies for arm(2), charge(2) and threshold(n_thrs)
        if self.do_range_full:
            if self.do_method_raw:
                self.full_raw_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

            if self.do_method_scaled:
                self.full_scaled_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

            if self.do_method_log:
                self.full_log_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

        if self.do_range_narrow:
            if self.do_method_raw:
                self.narrow_raw_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

            if self.do_method_scaled:
                self.narrow_scaled_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

            if self.do_method_log:
                self.narrow_log_instance = \
                    [[{k:0 for k in self.true_thrs} for j in range(2)]
                     for i in range(2)]

        self.status = 0
    # end def __init__()


    def get_data(self):
        """read data file and fill gpr_instance objects"""

        def normalize_vs_wness(array,uncert):
            '''helper function to normalize entries within each dw23 slice

            should result in a flat distribution with respect to wness
            '''
            sum_dict = dict()
            for wness in np.unique(array[:, 1]):
                sum_for_slices = np.sum(
                    array[array[:, 1]==np.around(wness,4)][:,2], axis=0)
                sum_dict[np.around(wness,4)] = sum_for_slices

            for i in range(0, len(array)):
                array[i, 2] /= sum_dict[np.around(array[i, 1], 4)]
                uncert[i] /= sum_dict[np.around(array[i, 1], 4)]
            return array[:, 2].reshape(len(array), 1), uncert, sum_dict
        # end def normalize_vs_wness()

        # file contains a list of entries describing bins from
        # a 2d dw23 vs wness yield histogram
        data = pd.read_csv(self.input_file,
                           sep=' ',
                           index_col=['index'],
                           usecols=['index', 'arm', 'charge',
                                    'wness_bin_center', 'dw23_bin_center',
                                    'entries'])

        # fill gpr_instance objects seaparately for each permutation of arm/charg/thr/method
        for arm in range(2):
            for charge in range(2):
                for threshold in self.true_thrs:
                    if self.do_range_full:
                        full_array = data[data['arm']==arm]\
                                    [data['charge']==charge]\
                                    [data['wness_bin_center']<threshold]\
                                    [data['wness_bin_center']>.1]\
                                    .as_matrix(['dw23_bin_center',
                                                'wness_bin_center',
                                                'entries'])

                        full_coordinates = \
                            full_array[:, [0, 1]].reshape(len(full_array),2)
                        full_values = \
                            full_array[:, 2].reshape(len(full_array), 1)
                        full_uncert = np.sqrt(full_values)

                    if self.do_range_narrow:
                        narrow_array = data[data['arm']==arm]\
                                [data['charge']==charge]\
                                [data['wness_bin_center']<threshold]\
                                [data['wness_bin_center']>.1]\
                                [data['dw23_bin_center']<.12]\
                                [data['dw23_bin_center']>-.12]\
                                .as_matrix(['dw23_bin_center',
                                            'wness_bin_center',
                                            'entries'])

                        narrow_coordinates = narrow_array[:, [0, 1]].reshape(
                            len(narrow_array), 2)
                        narrow_values = \
                            narrow_array[:, 2].reshape(len(narrow_array), 1)
                        narrow_uncert = np.sqrt(narrow_values)

                    if self.do_method_raw:
                        if self.do_range_full:
                            self.full_raw_instance[arm][charge][threshold] = \
                                gpr_instance('raw',
                                             'full',
                                             full_coordinates,
                                             full_values,
                                             full_uncert)

                        if self.do_range_narrow:
                            self.narrow_raw_instance[arm][charge][threshold] = \
                                gpr_instance('raw',
                                             'narrow',
                                             narrow_coordinates,
                                             narrow_values,
                                             narrow_uncert)

                    if self.do_method_scaled:
                        if self.do_range_full:
                            scaled_full_values, \
                                scaled_full_uncert, \
                                full_scale_factors = normalize_vs_wness(
                                    full_array, full_uncert)
                            self.full_scaled_instance[arm][charge][threshold] =\
                                    gpr_instance('scaled',
                                                 'full',
                                                 full_coordinates,
                                                 scaled_full_values,
                                                 scaled_full_uncert)
                            self.full_scaled_instance[arm][charge][threshold].\
                                store_scale_factors(full_scale_factors)
                            self.write_scale_factors(arm,
                                                     charge,
                                                     threshold,
                                                     'full',
                                                     full_scale_factors,
                                                     full_array)

                        if self.do_range_narrow:
                            scaled_narrow_values, \
                                scaled_narrow_uncert, \
                                narrow_scale_factors = normalize_vs_wness(
                                    narrow_array, narrow_uncert)
                            self.narrow_scaled_instance[arm][charge][threshold]\
                                = gpr_instance('scaled',
                                               'narrow',
                                               narrow_coordinates,
                                               scaled_narrow_values,
                                               scaled_narrow_uncert)
                            self.narrow_scaled_instance[arm][charge][threshold]\
                                .store_scale_factors(narrow_scale_factors)
                            self.write_scale_factors(arm,
                                                     charge,
                                                     threshold,
                                                     'narrow',
                                                     narrow_scale_factors,
                                                     narrow_array)



                    if self.do_method_log:
                        if self.do_range_full:
                            # log_full_values = np.log1p(full_values)
                            # log_full_values[np.invert(
                            #    np.isfinite(full_values))]=0.

                            # log_full_uncert = np.log(full_uncert)
                            # log_full_uncert[np.invert(
                            #    np.isfinite(full_uncert))]=0.

                            log_full_values = np.sqrt(full_values)
                            log_full_uncert = np.sqrt(full_uncert)

                            self.full_log_instance[arm][charge][threshold] = \
                                    gpr_instance('log',
                                                 'full',
                                                 full_coordinates,
                                                 log_full_values,
                                                 log_full_uncert)

                        if self.do_range_narrow:
                            #log_narrow_values = np.log1p(narrow_values)
                            #log_narrow_values[np.invert(
                            #    np.isfinite(narrow_values))]=0.

                            #log_narrow_uncert = np.log(narrow_uncert)
                            #log_narrow_uncert[np.invert(
                            #    np.isfinite(narrow_uncert))]=0.

                            log_narrow_values = np.sqrt(narrow_values)
                            log_narrow_uncert = np.sqrt(narrow_uncert)

                            self.narrow_log_instance[arm][charge][threshold] = \
                                    gpr_instance('log',
                                                 'narrow',
                                                 narrow_coordinates,
                                                 log_narrow_values,
                                                 log_narrow_uncert)

        # end gpr_instance filling permutation loop
    # end def get_data()

    def get_dw23_wness_coords(self,
                              dobincenter,
                              dw23bins,
                              dw23lower,
                              dw23upper,
                              wbins,
                              wlower,wupper):
        '''
        first argument is int(bool) on whether to
        takes integer/float input on spacing/limits for the data point grid
        outputs a 2d array that is a grid of x values for making predictions
        '''

        dw23increment = (dw23upper - dw23lower) / dw23bins
        wincrement = (wupper - wlower) / wbins

        if(dobincenter):
            dw23lower = dw23lower + dw23increment / 2
            dw23upper = dw23upper + dw23increment / 2
            wlower = wlower + wincrement / 2
            wupper = wupper + wincrement / 2

        dw23range = np.arange(dw23lower, dw23upper, dw23increment)
        wrange = np.arange(wlower, wupper, wincrement)

        first_time = 1
        for w in wrange:
            for dw23 in dw23range:
                if first_time:
                    predict_coords = np.array([[dw23, w]])
                    first_time = 0
                else:
                    predict_coords = np.concatenate((predict_coords,
                                                     [[dw23,w]]),
                                                    axis=0)

        return predict_coords

    def run_gpr(self):
        predict_coords = self.get_dw23_wness_coords(1, 120., -.3, .3,
                                                    200., 0., 1.)

        for threshold in self.true_thrs:
            for arm in range(1, 2):#changed from (2) to (1)
                for charge in range(2):
                    if self.do_range_full:
                        if self.do_method_raw:
                            print 'a%dc%dt%.1f raw create models...' % \
                                (arm, charge, threshold)
                            self.full_raw_instance[arm][charge][threshold].\
                                create_kernel_model()
                            print 'a%dc%dt%.1f raw optimize models...' % \
                                (arm, charge, threshold)
                            self.full_raw_instance[arm][charge][threshold].\
                                optimize_model()
                            print 'a%dc%dt%.1f raw generate predictions...' % \
                                (arm, charge, threshold)
                            self.full_raw_instance[arm][charge][threshold].\
                                generate_predictions(predict_coords)
                            print 'a%dc%dt%.1f raw write results...' % \
                                (arm, charge, threshold)
                            self.full_raw_instance[arm][charge][threshold].\
                                write_results(arm, charge, threshold)
                            print 'a%dc%dt%.1f raw clear memory...' % \
                                (arm, charge, threshold)
                            self.full_raw_instance[arm][charge][threshold].\
                                clear_memory()
                            del self.full_raw_instance[arm][charge][threshold]
                        if self.do_method_scaled:
                            print 'a%dc%dt%.1f scaled create models...' % \
                                (arm, charge, threshold)
                            self.full_scaled_instance[arm][charge][threshold].\
                                create_kernel_model()
                            print 'a%dc%dt%.1f scaled optimize models...' % \
                                (arm, charge, threshold)
                            self.full_scaled_instance[arm][charge][threshold].\
                                optimize_model()
                            print 'a%dc%dt%.1f scaled generate predictions...' \
                                % (arm, charge, threshold)
                            self.full_scaled_instance[arm][charge][threshold].\
                                generate_predictions(predict_coords)
                            self.full_scaled_instance[arm][charge][threshold].\
                                apply_rescaling()
                            print 'a%dc%dt%.1f scaled write results...' % \
                                (arm, charge, threshold)
                            self.full_scaled_instance[arm][charge][threshold].\
                                write_results(arm, charge, threshold)
                            print 'a%dc%dt%.1f scaled clear memory...' % \
                                (arm, charge, threshold)
                            self.full_scaled_instance[arm][charge][threshold].\
                                clear_memory()
                            del self.full_scaled_instance[arm]\
                                [charge][threshold]
                        if self.do_method_log:
                            print 'a%dc%dt%.1f log create models...' % \
                                (arm, charge, threshold)
                            self.full_log_instance[arm][charge][threshold].\
                                create_kernel_model()
                            print 'a%dc%dt%.1f log optimize models...' % \
                                (arm,charge,threshold)
                            self.full_log_instance[arm][charge][threshold].\
                                optimize_model()
                            print 'a%dc%dt%.1f log generate predictions...' % \
                                (arm, charge, threshold)
                            self.full_log_instance[arm][charge][threshold].\
                                generate_predictions(predict_coords)
                            self.full_log_instance[arm][charge][threshold].\
                                apply_rescaling()
                            print 'a%dc%dt%.1f log write results...' % \
                                (arm, charge, threshold)
                            self.full_log_instance[arm][charge][threshold].\
                                write_results(arm, charge, threshold)
                            print 'a%dc%dt%.1f log clear memory...' % \
                                (arm, charge, threshold)
                            self.full_log_instance[arm][charge][threshold].\
                                clear_memory()
                            del self.full_log_instance[arm][charge][threshold]
                    if self.do_range_narrow:
                        if self.do_method_raw:
                            print 'a%dc%dt%.1f raw create models...' % \
                                (arm, charge, threshold)
                            self.narrow_raw_instance[arm][charge][threshold].\
                                create_kernel_model()
                            print 'a%dc%dt%.1f raw optimize models...'%(arm,charge,threshold)
                            self.narrow_raw_instance[arm][charge][threshold].optimize_model()
                            print 'a%dc%dt%.1f raw generate predictions...'%(arm,charge,threshold)
                            self.narrow_raw_instance[arm][charge][threshold].generate_predictions(predict_coords)
                            print 'a%dc%dt%.1f raw write results...'%(arm,charge,threshold)
                            self.narrow_raw_instance[arm][charge][threshold].write_results(arm,charge,threshold)
                            print 'a%dc%dt%.1f raw clear memory...'%(arm,charge,threshold)
                            self.narrow_raw_instance[arm][charge][threshold].clear_memory()
                            del self.narrow_raw_instance[arm][charge][threshold]
                        if self.do_method_scaled:
                            print 'a%dc%dt%.1f scaled create models...'%(arm,charge,threshold)
                            self.narrow_scaled_instance[arm][charge][threshold].create_kernel_model()
                            print 'a%dc%dt%.1f scaled optimize models...'%(arm,charge,threshold)
                            self.narrow_scaled_instance[arm][charge][threshold].optimize_model()
                            print 'a%dc%dt%.1f scaled generate predictions...'%(arm,charge,threshold)
                            self.narrow_scaled_instance[arm][charge][threshold].generate_predictions(predict_coords)
                            self.narrow_scaled_instance[arm][charge][threshold].apply_rescaling()
                            print 'a%dc%dt%.1f scaled write results...'%(arm,charge,threshold)
                            self.narrow_scaled_instance[arm][charge][threshold].write_results(arm,charge,threshold)
                            print 'a%dc%dt%.1f scaled clear memory...'%(arm,charge,threshold)
                            self.narrow_scaled_instance[arm][charge][threshold].clear_memory()
                            del self.narrow_scaled_instance[arm][charge][threshold]
                        if self.do_method_log:
                            print 'a%dc%dt%.1f log create models...'%(arm,charge,threshold)
                            self.narrow_log_instance[arm][charge][threshold].create_kernel_model()
                            print 'a%dc%dt%.1f log optimize models...'%(arm,charge,threshold)
                            self.narrow_log_instance[arm][charge][threshold].optimize_model()
                            print 'a%dc%dt%.1f log generate predictions...'%(arm,charge,threshold)
                            self.narrow_log_instance[arm][charge][threshold].generate_predictions(predict_coords)
                            self.narrow_log_instance[arm][charge][threshold].apply_rescaling()
                            print 'a%dc%dt%.1f log write results...'%(arm,charge,threshold)
                            self.narrow_log_instance[arm][charge][threshold].write_results(arm,charge,threshold)
                            print 'a%dc%dt%.1f log clear memory...'%(arm,charge,threshold)
                            self.narrow_log_instance[arm][charge][threshold].clear_memory()
                            del self.narrow_log_instance[arm][charge][threshold]

        ####

    def write_scale_factors(self,arm,charge,threshold,dw23range,scale_factors,array):
        filename = 'output/checks/scale_factors_input_thr%.1f_%s_a%d_c%d.txt'%(threshold,dw23range,arm,charge)
        print('writing scale factors to file: %s'%filename)
        fout = open(filename, 'w')
        fout.write("index scale_factor scale_count dw23_bin_center wness_bin_center entries\n")
        for i in range(0,len(array)):
            fout.write("%5d %12f %11.0f %15f %16f %7f\n"%(i,
                        1/scale_factors[np.around(array[i,1],4)],
                        scale_factors[np.around(array[i,1],4)],
                        array[i,0],
                        array[i,1],
                        array[i,2]))

    def testprint(self):
        for arm in range(2):
            for charge in range(2):
                for threshold in self.true_thrs:
                    print 'arm %d charge %d thr %f'%(arm,charge,threshold)
                    if self.do_range_full:
                        if self.do_method_raw:
                            self.full_raw_instance[arm][charge][threshold].testprint()
                        if self.do_method_scaled:
                            self.full_scaled_instance[arm][charge][threshold].testprint()
                        if self.do_method_log:
                            self.full_log_instance[arm][charge][threshold].testprint()
                    if self.do_range_narrow:
                        if self.do_method_raw:
                            self.narrow_raw_instance[arm][charge][threshold].testprint()
                        if self.do_method_scaled:
                            self.narrow_scaled_instance[arm][charge][threshold].testprint()
                        if self.do_method_log:
                            self.narrow_log_instance[arm][charge][threshold].testprint()

#end gpr_obj class definition

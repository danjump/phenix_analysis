import numpy as np
import pandas as pd

from gpr_instance import gpr_instance

class gpr_analysis_obj:
    """a master class for the grp analysis
    
    contains functions and variables to facilitate 
    the flow of executing the analysisi
    """

    def __init__(self,args_dict):
        """initialize class object with analysis options"""
        self.do_method_raw =    args_dict['do_method_raw'   ]
        self.do_method_scaled = args_dict['do_method_scaled']
        self.do_method_log =    args_dict['do_method_log'   ]
        self.do_threshold =     args_dict['do_threshold'    ]
        self.input_file =       args_dict['input_file'      ]

        # generate array of thresholds that we want to do (True ones only)
        self.true_thrs = self.do_threshold[self.do_threshold[:,1]==True][:,0]

        # zero-initialize arrays that will hold gpr_instance objects
        # 3 indecies for arm(2), charge(2) and threshold(n_thrs)
        if self.do_method_raw:
            self.raw_instance = [[{k:0 for k in self.true_thrs} for j in range(2)] for i in range(2)]
        
        if self.do_method_scaled:
            self.scaled_instance = [[{k:0 for k in self.true_thrs} for j in range(2)] for i in range(2)]
        
        if self.do_method_log:
            self.log_instance = [[{k:0 for k in self.true_thrs} for j in range(2)] for i in range(2)]
        
        self.status = 0
    # end def __init__()

    def get_data(self):
        """read data file and fill gpr_instance objects"""
        
        def normalize_vs_wness(array,uncert):
            '''helper function to normalize entries within each dw23 slice

            should result in a flat distribution with respect to wness
            '''
            sum_dict = dict()
            for wness in np.unique(array[:,1]):
                sum_for_slices = np.sum(array[array[:,1]==np.around(wness,2)][:,2],axis=0)
                sum_dict[np.around(wness,2)] = sum_for_slices
                
            for i in range(0,len(array)):
                array[i,2] /= sum_dict[np.around(array[i,1],2)]
                uncert[i] /= sum_dict[np.around(array[i,1],2)]
            return array[:,2].reshape(len(array),1), uncert, sum_dict
        # end def normalize_vs_wness()

        # file contains a list of entries describing bins from
        # a 2d dw23 vs wness yield histogram
        data = pd.read_csv(self.input_file,
                        sep=' ',
                        index_col=['index'],
                        usecols=['index','arm','charge','wness_bin_center',
                                 'dw23_bin_center','entries'])

        # fill gpr_instance objects seaparately for each permutation of arm/charg/thr/method
        for arm in range(2):
            for charge in range(2):
                for threshold in self.true_thrs:
                    full_array = data[data['arm']==arm]\
                                [data['charge']==charge]\
                                [data['wness_bin_center']<threshold]\
                                [data['wness_bin_center']>.1]\
                                .as_matrix(['dw23_bin_center','wness_bin_center','entries'])

                    coordinates = full_array[:,[0,1]].reshape(len(full_array),2) 
                    values = full_array[:,2].reshape(len(full_array),1)
                    uncert = np.sqrt(values)

                    if self.do_method_raw:
                        self.raw_instance[arm][charge][threshold] = gpr_instance('raw',coordinates,values,uncert)

                    if self.do_method_scaled:
                        scaled_values,scaled_uncert,scale_factors = normalize_vs_wness(full_array,uncert)
                        self.scaled_instance[arm][charge][threshold] = \
                                gpr_instance('scaled',coordinates,scaled_values,scaled_uncert)
                        self.scaled_instance[arm][charge][threshold].store_scale_factors(scale_factors)

                    if self.do_method_log:
                        log_values = np.log(values+2)
                        log_values[np.invert(np.isfinite(values))]=0.

                        log_uncert = np.log(uncert)
                        log_uncert[np.invert(np.isfinite(uncert))]=0.

                        self.log_instance[arm][charge][threshold] = \
                                gpr_instance('log',coordinates,log_values,log_uncert)
        # end gpr_instance filling permutation loop 
    # end def get_data()
    
    def get_dw23_wness_coords(self,dobincenter,dw23bins,dw23lower,dw23upper,wbins,wlower,wupper):
        '''
        first argument is int(bool) on whether to 
        takes integer/float input on spacing/limits for the data point grid
        outputs a 2d array that is a grid of x values for making predictions
        '''

        dw23increment = (dw23upper-dw23lower)/dw23bins
        wincrement = (wupper-wlower)/wbins

        if(dobincenter):
            dw23lower = dw23lower+dw23increment/2
            dw23upper = dw23upper+dw23increment/2
            wlower = wlower+wincrement/2
            wupper = wupper+wincrement/2
        
        dw23range = np.arange(dw23lower,dw23upper,dw23increment)
        wrange = np.arange(wlower,wupper,wincrement)

        first_time = 1
        for dw23 in dw23range:
            for w in wrange:
                if first_time:
                    predict_coords = np.array([[dw23,w]])
                    first_time = 0
                else:
                    predict_coords = np.concatenate((predict_coords,[[dw23,w]]),axis=0)
                    
        return predict_coords
    
    def run_gpr(self):
        predict_coords = self.get_dw23_wness_coords(1, 120.,-.3,.3, 200.,0.,1.)
        
        for threshold in self.true_thrs:
            for arm in range(2):
                for charge in range(2):
                    if self.do_method_raw:
                        print 'a%dc%dt%.1f raw create models...'%(arm,charge,threshold)
                        self.raw_instance[arm][charge][threshold].create_kernel_model()
                        print 'a%dc%dt%.1f raw optimize models...'%(arm,charge,threshold)
                        self.raw_instance[arm][charge][threshold].optimize_model()
                        print 'a%dc%dt%.1f raw generate predictions...'%(arm,charge,threshold)
                        self.raw_instance[arm][charge][threshold].generate_predictions(predict_coords)
                        print 'a%dc%dt%.1f raw write results...'%(arm,charge,threshold)
                        results = self.raw_instance[arm][charge][threshold].get_results()
                        self.write_results(results,predict_coords,'raw',arm,charge,threshold)
                        print 'a%dc%dt%.1f raw clear memory...'%(arm,charge,threshold)
                        self.raw_instance[arm][charge][threshold].clear_memory()
                        del self.raw_instance[arm][charge][threshold]
                    if self.do_method_scaled:
                        print 'a%dc%dt%.1f scaled create models...'%(arm,charge,threshold)
                        self.scaled_instance[arm][charge][threshold].create_kernel_model()
                        print 'a%dc%dt%.1f scaled optimize models...'%(arm,charge,threshold)
                        self.scaled_instance[arm][charge][threshold].optimize_model()
                        print 'a%dc%dt%.1f scaled generate predictions...'%(arm,charge,threshold)
                        self.scaled_instance[arm][charge][threshold].generate_predictions(predict_coords)
                        self.scaled_instance[arm][charge][threshold].apply_rescaling()
                        print 'a%dc%dt%.1f scaled write results...'%(arm,charge,threshold)
                        results = self.scaled_instance[arm][charge][threshold].get_results()
                        self.write_results(results,predict_coords,'scaled',arm,charge,threshold)
                        print 'a%dc%dt%.1f scaled clear memory...'%(arm,charge,threshold)
                        self.scaled_instance[arm][charge][threshold].clear_memory()
                        del self.scaled_instance[arm][charge][threshold]
                    if self.do_method_log:
                        print 'a%dc%dt%.1f log create models...'%(arm,charge,threshold)
                        self.log_instance[arm][charge][threshold].create_kernel_model()
                        print 'a%dc%dt%.1f log optimize models...'%(arm,charge,threshold)
                        self.log_instance[arm][charge][threshold].optimize_model()
                        print 'a%dc%dt%.1f log generate predictions...'%(arm,charge,threshold)
                        self.log_instance[arm][charge][threshold].generate_predictions(predict_coords)
                        print 'a%dc%dt%.1f log write results...'%(arm,charge,threshold)
                        results = self.log_instance[arm][charge][threshold].get_results()
                        self.write_results(results,predict_coords,'log',arm,charge,threshold)
                        print 'a%dc%dt%.1f log clear memory...'%(arm,charge,threshold)
                        self.log_instance[arm][charge][threshold].clear_memory()
                        del self.log_instance[arm][charge][threshold]
                        
        ####
    def write_results(self,results,predict_coords,method,arm,charge,threshold):
        mean =        results[0]
        var =         results[1]
        mean_uncert = results[2]
        var_uncert =  results[3]
        filename = 'output/gpr_predicted_points_%s_thr%.1f_a%d_c%d.txt'%(method,threshold,arm,charge)
        fout = open(filename, 'w')
        fout.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")
        for i in range(0,len(mean)):
            fout.write("%d %d %d %f %f %f %f %f %f %f\n"%(i,arm,charge,\
                np.around(threshold,1),\
                predict_coords[i][0],\
                predict_coords[i][1],\
                mean[i],\
                var[i],\
                mean_uncert[i],\
                var_uncert[i]))
        fout.close()
    # end def write_results

    def write_scale_factors(self,arm,charge,threshold,scale_factors,full_array):
        filename = 'output/dw23_slice_scale_factors_thr%.1f_a%d_c%d.txt'%(threshold,arm,charge)
        print('writing scale factors to file: %s'%filename)
        fout = open(filename, 'w')
        fout.write("index scale_factor scale_count dw23_bin_center wness_bin_center entries\n")
        for i in range(0,len(full_array)):
            fout.write("%5d %12f %11.0f %15f %16f %7f\n"%(i,
                        1/scale_factors[np.around(full_array[i,1],2)],
                        scale_factors[np.around(full_array[i,1],2)],
                        full_array[i,0],
                        full_array[i,1],
                        full_array[i,2]))

    def testprint(self):
        for arm in range(2):
            for charge in range(2):
                for threshold in self.true_thrs:
                    print 'arm %d charge %d thr %f'%(arm,charge,threshold)
                    if self.do_method_raw:
                        self.raw_instance[arm][charge][threshold].testprint()
                    if self.do_method_scaled:
                        self.scaled_instance[arm][charge][threshold].testprint()
                    if self.do_method_log:
                        self.log_instance[arm][charge][threshold].testprint()

#end gpr_obj class definition

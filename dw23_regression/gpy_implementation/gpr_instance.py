import GPy
import numpy as np


class gpr_instance:

    def __init__(self,method,dw23range,coords,values,uncert):
        self.coords = coords
        self.values = values
        self.uncert = uncert
        self.method = method
        self.dw23range = dw23range
        if method == 'raw':
            self.dw23variance = 100
            self.wvariance = 100
        elif method == 'scaled':
            self.dw23variance = 0.1
            self.wvariance = 0.01
        elif method == 'log':
            self.dw23variance = 10
            self.wvariance = 10


    def testprint(self):
        print self.method,self.dw23range,self.coords[15],self.values[15],self.uncert[15]

    def create_kernel_model(self):
        # create model without input uncertainty
        self.kerdw23 = GPy.kern.RBF(input_dim=1,variance=self.dw23variance,lengthscale=.05,active_dims=[0])
        self.kerw = GPy.kern.RBF(input_dim=1,variance=self.wvariance,lengthscale=.05,active_dims=[1])
        self.ker = self.kerdw23 * self.kerw

        self.m = GPy.models.GPRegression(self.coords,self.values,self.ker)

        # create model using input uncertainty:
        self.kerdw23uncert = GPy.kern.RBF(input_dim=1,variance=self.dw23variance,lengthscale=.05,active_dims=[0])
        self.kerwuncert = GPy.kern.RBF(input_dim=1,variance=self.wvariance,lengthscale=.05,active_dims=[1])
        self.kerfix = GPy.kern.Fixed(2,GPy.util.linalg.tdot(self.uncert))
        self.keruncert = self.kerdw23uncert * self.kerwuncert + self.kerfix

        self.m_uncert = GPy.models.GPRegression(self.coords,self.values,self.keruncert)
        self.m_uncert.Gaussian_noise.fix(1e-6)

    def optimize_model(self):
        self.m.optimize(max_f_eval = 1000)
        #print(self.m)
        self.m_uncert.optimize(max_f_eval = 1000)
        #print(self.m_uncert)

    def generate_predictions(self,predict_coords):
        self.predict_coords = predict_coords
        self.mean, self.var = self.m.predict(predict_coords)
        self.mean_uncert, self.var_uncert = self.m_uncert.predict(
            predict_coords,
            kern = self.kerdw23uncert.copy()*self.kerwuncert.copy())

    def write_results(self,arm,charge,threshold):
        #write final prediction results:
        filename = 'output/ralf_hadsim/gpr_results/gpr_predicted_points_%s_%s_thr%.1f_a%d_c%d.txt'%(self.method,self.dw23range,threshold,arm,charge)
        fout = open(filename, 'w')
        fout.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")
        for i in range(0,len(self.mean)):
            fout.write("%d %d %d %f %f %f %f %f %f %f\n"%(i,arm,charge,\
                np.around(threshold,1),\
                self.predict_coords[i][0],\
                self.predict_coords[i][1],\
                self.mean[i],\
                self.var[i],\
                self.mean_uncert[i],\
                self.var_uncert[i]))
        fout.close()

        #write applied_scale_factors and scaled predictions
        if self.method=='scaled':
            filename = 'output/ralf_hadsim/checks/scale_factors_output_%s_%s_thr%.1f_a%d_c%d.txt'%(self.method,self.dw23range,threshold,arm,charge)
            fout = open(filename,'w')
            fout.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")
            for i in range(0,len(self.mean)):
                fout.write("%d %d %d %f %f %f %f %f %f %f %f\n"%(i,arm,charge,\
                    np.around(threshold,1),\
                    self.predict_coords[i][0],\
                    self.predict_coords[i][1],\
                    self.mean[i],\
                    self.var[i],\
                    self.mean_uncert[i],\
                    self.var_uncert[i],\
                    self.applied_scale_factors[np.around(self.predict_coords[i][1],4)]))
            fout.close()
        #end if method=='scaled'

    # end def write_results

    def clear_memory(self):
        del self.kerdw23
        del self.kerw
        del self.ker
        del self.m
        del self.kerdw23uncert
        del self.kerwuncert
        del self.kerfix
        del self.keruncert
        del self.m_uncert
        del self.mean
        del self.var
        del self.mean_uncert
        del self.var_uncert
        del self.coords
        del self.values
        del self.uncert
        del self.predict_coords
        if self.method == 'scaled':
            del self.scale_factors
            del self.applied_scale_factors

    def get_results(self):
        return self.mean, self.var, self.mean_uncert, self.var_uncert

    def store_scale_factors(self,scale_factors_dict):
        self.scale_factors = scale_factors_dict

    def apply_rescaling(self):
        self.rescaled_flag = True
        if self.method == 'scaled':
            scale_factors_list = list()
            for i in range(0,len(self.predict_coords)):
                wness = self.predict_coords[i,1]
                scale_factor = self.calculate_scale_factor(wness)
                self.mean[i] *= scale_factor
                self.var[i] *= scale_factor
                self.mean_uncert[i] *= scale_factor
                self.var_uncert[i] *= scale_factor
                scale_factors_list.append([wness,scale_factor])

            scale_factor_array = np.array(scale_factors_list)
            dummy,index = np.unique(scale_factor_array[:,0],return_index=True)
            self.applied_scale_factors = dict()
            for i in index:
                wness = scale_factor_array[i,0]
                scale_factor = scale_factor_array[i,1]
                self.applied_scale_factors[np.around(wness,4)] = scale_factor
        if self.method == 'log':
            #for i in range(0,len(self.predict_coords)):
                #self.mean[i] = np.expm1(self.mean[i])
                #self.var[i] = np.exp(self.var[i])
                #self.mean_uncert[i] = np.expm1(self.mean_uncert[i])
                #self.var_uncert[i] = np.exp(self.var_uncert[i])
            self.mean = np.power(self.mean,2)
            self.var = np.power(self.var,2)
            self.mean_uncert = np.power(self.mean_uncert,2)
            self.var_uncert = np.power(self.var_uncert,2)

    #end def apply_rescaling()

    def calculate_scale_factor(self,wness):
        wness_scale_points = np.unique(self.coords[:,1])
        if wness <= wness_scale_points[0]:
            scale_factor = self.scale_factors[np.around(wness_scale_points[0],4)]
        elif wness >= wness_scale_points[-1]:
            scale_factor = self.scale_factors[np.around(wness_scale_points[-1],4)]
        else:
            for i in range(0,len(wness_scale_points)-1):
                if wness > wness_scale_points[i+1]:
                    continue
                else:
                    if np.isclose(wness,wness_scale_points[i]):
                        scale_factor = self.scale_factors[np.around(wness_scale_points[i],4)]
                    else:
                        lower_wness = wness_scale_points[i]
                        lower_factor = self.scale_factors[lower_wness]
                        upper_wness = wness_scale_points[i+1]
                        upper_factor = self.scale_factors[upper_wness]

                        lower_weight = (wness-lower_wness)/(upper_wness-lower_wness)
                        upper_weight = (upper_wness-wness)/(upper_wness-lower_wness)

                        scale_factor = lower_factor*lower_weight + upper_factor*upper_weight
                    break
                #end else
            #end for i ....
        #end else
        return scale_factor
    #end def calculate_scale_factor()








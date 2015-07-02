import GPy

class gpr_instance:
    
    def __init__(self,method,coords,values,uncert):
        self.coords = coords
        self.values = values
        self.uncert = uncert
        self.method = method
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
        print self.method,self.coords[15],self.values[15],self.uncert[15]

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
        self.mean, self.var = self.m.predict(predict_coords)
        self.mean_uncert, self.var_uncert = self.m_uncert.predict(
            predict_coords,
            kern = self.kerdw23uncert.copy()*self.kerwuncert.copy())
        
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
        
    def get_results(self):
        return self.mean, self.var, self.mean_uncert, self.var_uncert 

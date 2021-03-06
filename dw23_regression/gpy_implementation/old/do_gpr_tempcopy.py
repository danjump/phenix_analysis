import pylab as pb
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import GPy as GPy
pb.ion()

#set running configuration variables:
do_method_raw=False
do_method_scaled=True
do_method_log=False

threshold = 2

# #------------------------------------
# # Read/Prepare input data
# #------------------------------------

# read from file
# file contains a list of entries describing bins from
# a 2d dw23 vs wness yield histogram
data = pd.read_csv('dw23_vs_wness_list.txt',
                sep=' ',
                index_col=['index'],
                usecols=['index','arm','charge','wness_bin_center',
                         'dw23_bin_center','entries'])

# select the arm/charge/wness region for the data points
def get_data_array(data,which,arm,charge,wthreshold):
    ''' simple function to filter events from dataframe and get array
    which=0 -> coordinates 2 column array
    which=1 -> values 1 column array
    which=2 -> coordinates and values 3 column array
    '''
    
    if which==1:
        cols = ['entries']
    elif which==0:
        cols = ['dw23_bin_center','wness_bin_center']
    elif which==2:
        cols = ['dw23_bin_center','wness_bin_center','entries']
        
    array = data[data['arm']==arm][data['charge']==charge][data['wness_bin_center']<wthreshold]            [data['wness_bin_center']>.1].as_matrix(cols)
        
    return array

#zero-initialize variables as 2,2,3 arrays for arm,charge,threshold
coordinates = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
values = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
full_array = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
scaled_values = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
scaled_uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

def normalize_vs_wness(array,uncert):
    '''
    normalizes entries within each dw23 slice.
    should result in a flat distribution with respect to wness
    '''
    sum_dict = dict()
    for dw23 in np.arange(-.29,.31,.02):
        dw23 = np.around(dw23,2)
        summed_vs_wness = np.sum(array[array[:,0]==dw23][:,2],axis=0)
        sum_dict[dw23] = summed_vs_wness
        
    for i in range(0,len(array)):
        array[i,2] /= sum_dict[array[i,0]]
        uncert[i] /= sum_dict[array[i,0]]
    return array[:,2].reshape(len(array),1), uncert

for arm in range(2):
    for charge in range(2):
        #assign un-modified data content
        coordinates[arm][charge][threshold] = get_data_array(data,0,arm,charge,.7+threshold/10.)
        values[arm][charge][threshold] = get_data_array(data,1,arm,charge,.7+threshold/10.)
        uncert[arm][charge][threshold] = np.sqrt(values[arm][charge][threshold])
        
        if(do_method_scaled):
            full_array[arm][charge][threshold] = get_data_array(data,2,arm,charge,.7+threshold/10.)
            scaled_values[arm][charge][threshold],scaled_uncert[arm][charge][threshold] = normalize_vs_wness(full_array[arm][charge][threshold],uncert[arm][charge][threshold])
            


# take the log of all values to make them compatible with the GPR fit
# note: add 2 to transform
if(do_method_log):
    log_values = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    log_uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    for arm in range(2):
        for charge in range(2):
            log_values[arm][charge][threshold] = np.log(values[arm][charge][threshold]+2)
            log_values[arm][charge][threshold][np.invert(np.isfinite(values[arm][charge][threshold]))]=0
            log_uncert[arm][charge][threshold] = np.log(uncert[arm][charge][threshold])
            log_uncert[arm][charge][threshold][np.invert(np.isfinite(uncert[arm][charge][threshold]))]=0.


# #-----------------------------------------
# # Do Gaussian Proccess Regression on Data
# #-----------------------------------------

# define kernels
if do_method_raw:    
    kerdw23 = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerw = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    ker = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    kerdw23uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerwuncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerfix = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    keruncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    m = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    m_uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    for arm in range(2):
        for charge in range(2):
            # create simple model:
            kerdw23[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[0])
            kerw[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[1])
            ker[arm][charge][threshold] = kerdw23[arm][charge][threshold] * kerw[arm][charge][threshold]
            
            m[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],values[arm][charge][threshold],ker[arm][charge][threshold])
            
            # create model using input uncertainty:
            kerdw23uncert[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[0])
            kerwuncert[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[1])
            kerfix[arm][charge][threshold] = GPy.kern.Fixed(2,GPy.util.linalg.tdot(uncert[arm][charge][threshold]))
            keruncert[arm][charge][threshold] = kerdw23uncert[arm][charge][threshold] * kerwuncert[arm][charge][threshold] + kerfix[arm][charge][threshold]
            
            m_uncert[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],values[arm][charge][threshold],keruncert[arm][charge][threshold])
            m_uncert[arm][charge][threshold].Gaussian_noise.fix(1e-6)
#end if do_method_raw

if do_method_scaled:
    kerdw23_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerw_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    ker_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    kerdw23uncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerwuncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerfix_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    keruncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    m_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    m_uncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    for arm in range(2):
        for charge in range(2):
            # create simple model:
            kerdw23_scaled[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=.1,lengthscale=.05,active_dims=[0])
            kerw_scaled[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=.01,lengthscale=.05,active_dims=[1])
            ker_scaled[arm][charge][threshold] = kerdw23_scaled[arm][charge][threshold] * kerw_scaled[arm][charge][threshold]
             
            m_scaled[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],scaled_values[arm][charge][threshold],ker_scaled[arm][charge][threshold])
            
            # create model using input uncertainty:
            kerdw23uncert_scaled[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=.1,lengthscale=.05,active_dims=[0])
            kerwuncert_scaled[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=.01,lengthscale=.05,active_dims=[1])
            kerfix_scaled[arm][charge][threshold] = GPy.kern.Fixed(2,GPy.util.linalg.tdot(scaled_uncert[arm][charge][threshold]))
            keruncert_scaled[arm][charge][threshold] = kerdw23uncert_scaled[arm][charge][threshold] * kerwuncert_scaled[arm][charge][threshold] + kerfix_scaled[arm][charge][threshold]
            
            m_uncert_scaled[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],scaled_values[arm][charge][threshold],keruncert_scaled[arm][charge][threshold])
            m_uncert_scaled[arm][charge][threshold].Gaussian_noise.fix(1e-6)
#end if do_method_scaled

if do_method_log:
    kerdw23_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerw_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    ker_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    kerdw23uncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerwuncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    kerfix_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    keruncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    m_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    m_uncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]


    for arm in range(2):
        for charge in range(2):
            # create simple model:
            kerdw23_log[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[0])
            kerw_log[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[1])
            ker_log[arm][charge][threshold] = kerdw23_log[arm][charge][threshold] * kerw_log[arm][charge][threshold]
            
            m_log[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],log_values[arm][charge][threshold],ker_log[arm][charge][threshold])
            
            # create model using input uncertainty:
            kerdw23uncert_log[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[0])
            kerwuncert_log[arm][charge][threshold] = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[1])
            kerfix_log[arm][charge][threshold] = GPy.kern.Fixed(2,GPy.util.linalg.tdot(log_uncert[arm][charge][threshold]))
            keruncert_log[arm][charge][threshold] = kerdw23uncert_log[arm][charge][threshold] * kerwuncert_log[arm][charge][threshold] + kerfix_log[arm][charge][threshold]
            
            m_uncert_log[arm][charge][threshold] = GPy.models.GPRegression(coordinates[arm][charge][threshold],log_values[arm][charge][threshold],keruncert_log[arm][charge][threshold])
            m_uncert_log[arm][charge][threshold].Gaussian_noise.fix(1e-6)
#end if do_method_log


# optimize kernel parameters (i.e. - "run" the GPR)
if do_method_raw:
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            print('Optimizing Raw Model...')
            m[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m[arm][charge][threshold])
            print('Optimizing Raw Model with uncertainty...')
            m_uncert[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m_uncert[arm][charge][threshold])
# end if do_method_raw

# optimize kernel parameters (i.e. - "run" the GPR)
if do_method_scaled:
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            print('Optimizing Scaled Model...')
            m_scaled[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m_scaled[arm][charge][threshold])
            print('Optimizing Scaled Model with uncertainty...')
            m_uncert_scaled[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m_uncert_scaled[arm][charge][threshold])
#end if do_method_scaled

# optimize kernel parameters (i.e. - "run" the GPR)
if do_method_log:
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            print('Optimizing Log Model...')
            m_log[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m_log[arm][charge][threshold])
            print('Optimizing Log Model with uncertainty...')
            m_uncert_log[arm][charge][threshold].optimize(max_f_eval = 1000)
            print(m_uncert_log[arm][charge][threshold])
#end if do_method_log

# #--------------------------------
# # Get GPR Results & Write to file
# #--------------------------------

# create a grid of dw23 & wness points to use
# for prediction points from our model
def get_dw23_wness_coords(dobincenter,dw23bins,dw23lower,dw23upper,wbins,wlower,wupper):
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

predict_coords = get_dw23_wness_coords(1, 120.,-.3,.3, 200.,0.,1.)

if do_method_raw:
    mean = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    mean_uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var_uncert = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    print 'Calculating predictions for raw models...'
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            #extract prediced mean/variance values from model
            mean[arm][charge][threshold], var[arm][charge][threshold] = m[arm][charge][threshold].predict(predict_coords)
            mean_uncert[arm][charge][threshold], var_uncert[arm][charge][threshold] =             m_uncert[arm][charge][threshold].predict(predict_coords,kernel=kerdw23uncert[arm][charge][threshold].copy()*kerwuncert[arm][charge][threshold].copy())
#end if do_method_raw

if do_method_scaled:
    mean_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    mean_uncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var_uncert_scaled = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    print 'Calculating predictions for scaled models...'
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            #extract prediced mean/variance values from model
            mean_scaled[arm][charge][threshold], var_scaled[arm][charge][threshold] = m_scaled[arm][charge][threshold].predict(predict_coords)
            mean_uncert_scaled[arm][charge][threshold], var_uncert_scaled[arm][charge][threshold] = m_uncert_scaled[arm][charge][threshold].predict(predict_coords)
#end if do_method_scaled

if do_method_log:
    mean_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    mean_uncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]
    var_uncert_log = [[[0 for x in range(3)] for x in range(2)] for x in range(2)]

    print 'Calculating predictions for log models...'
    for arm in range(2):
        for charge in range(2):
            print '\n\nArm%d Charge%d Threshold%.1f\n'%(arm,charge,.7+threshold/10.)
            #extract prediced mean/variance values from model
            mean_log[arm][charge][threshold], var_log[arm][charge][threshold] = m_log[arm][charge][threshold].predict(predict_coords)
            mean_uncert_log[arm][charge][threshold], var_uncert_log[arm][charge][threshold] = m_uncert_log[arm][charge][threshold].predict(predict_coords)

            # exectue exponential function on results to un-do the log() we did before fitting
            mean_log[arm][charge][threshold]=np.exp(mean_log[arm][charge][threshold])-2
            var_log[arm][charge][threshold]=np.exp(var_log[arm][charge][threshold])
            mean_uncert_log[arm][charge][threshold]=np.exp(mean_uncert_log[arm][charge][threshold])-2
            var_uncert_log[arm][charge][threshold]=np.exp(var_uncert_log[arm][charge][threshold])


#write results to file
if do_method_raw:
    fout_raw = open('output/predicted_dw23_wness_points_raw.txt', 'w')
    fout_raw.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")

if do_method_scaled:
    fout_scaled = open('output/predicted_dw23_wness_points_scaled.txt', 'w')
    fout_scaled.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")

if do_method_log:
    fout_scaled = open('output/predicted_dw23_wness_points_log.txt', 'w')
    fout_scaled.write("index arm charge threshold dw23_bin_center wness_bin_center nouncert_value nouncert_uncertainty withuncert_value withuncert_uncertainty\n")

for arm in range(2):
    for charge in range(2):
        if do_method_raw:
            for i in range(0,len(mean[arm][charge][threshold])):
                fout_raw.write("%d %d %d %f %f %f %f %f\n"%(i,arm,charge,np.around(.7+threshold/10.,1),predict_coords[i][0],predict_coords[i][1],mean[arm][charge][threshold][i],var[arm][charge][threshold][i],mean_uncert[arm][charge][threshold][i],var_uncert[arm][charge][threshold][i]))
        
        if do_method_scaled:
            for i in range(0,len(mean[arm][charge][threshold])):
                fout_scaled.write("%d %d %d %f %f %f %f %f\n"%(i,arm,charge,np.around(.7+threshold/10.,1),predict_coords[i][0],predict_coords[i][1],mean[arm][charge][threshold][i],var[arm][charge][threshold][i],mean_uncert[arm][charge][threshold][i],var_uncert[arm][charge][threshold][i]))

        if do_method_log:
            for i in range(0,len(mean[arm][charge][threshold])):
                fout_log.write("%d %d %d %f %f %f %f %f\n"%(i,arm,charge,np.around(.7+threshold/10.,1),predict_coords[i][0],predict_coords[i][1],mean[arm][charge][threshold][i],var[arm][charge][threshold][i],mean_uncert[arm][charge][threshold][i],var_uncert[arm][charge][threshold][i]))
            


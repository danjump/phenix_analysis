
# coding: utf-8

# In[64]:

import pylab as pb
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import GPy
pb.ion()


# In[3]:




# read dw23/wness entries
data = pd.read_csv('dw23_vs_wness_list.txt',
                sep=' ',
                index_col=['index'],
                usecols=['index','arm','charge','wness_bin_center','dw23_bin_center','entries'])


coordinates = data[data['arm']==0][data['charge']==0][data['wness_bin_center']<.9][data['wness_bin_center']>.1].as_matrix(['dw23_bin_center','wness_bin_center'])

values = data[data['arm']==0][data['charge']==0][data['wness_bin_center']<.9][data['wness_bin_center']>.1].as_matrix(['entries'])


# In[4]:

# define kernel
ker1 = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[0])
ker2 = GPy.kern.RBF(input_dim=1,variance=100,lengthscale=.05,active_dims=[1])

ker = ker1 * ker2

# create simple GP model
m = GPy.models.GPRegression(coordinates,values,ker)

# optimize and plot
m.optimize(max_f_eval = 1000)
print(m)


# In[115]:

test = 1

dw23range = np.arange(-.3,.3,.02)
wrange = np.arange(0.1,.9,.02)

for dw23 in dw23range:
    for w in wrange:
        if test:
            predicts = np.array([[dw23,w]])
            test = 0
        else:
            predicts = np.concatenate((predicts,[[dw23,w]]),axis=0)
            
wrange,dw23range = np.meshgrid(wrange,dw23range)
        
        
mean, var = m.predict(predicts)
mean = mean.reshape(dw23range.shape)
var = var.reshape(dw23range.shape)

mean.shape


# In[122]:

scatterx,scattery = np.hsplit(coordinates,2)
scatterx = scatterx.reshape(40,30).T
scattery = scattery.reshape(40,30).T
scatterz = values.reshape(40,30).T


# In[112]:


m.plot()

'''
figure, axes = plt.subplots(4,1)
for ax,y in zip(axes,[.2,.4,.6,.8]):
  m.plot(fixed_inputs=[(1,y)],ax=ax, which_data_rows=slice(0,30))
  plt.xlabel('wness prediction %i'%y)'''
upperBand = mean + var
lowerBand = mean - var

fig=plt.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(dw23range,wrange,mean,colors='blue')
ax.plot_wireframe(dw23range,wrange,upperBand,colors='red')
ax.plot_wireframe(dw23range,wrange,lowerBand,colors='red')
ax.scatter(scatterx,scattery,scatterz)
ax.set_xlabel('dw23')
ax.set_ylabel('wness')
plt.show()


# In[124]:

fig=plt.figure()
ax = fig.gca(projection='3d')
#ax.plot_wireframe(dw23range,wrange,mean-scatterz,colors='blue')
ax.scatter(dw23range,wrange,mean-scatterz)
ax.plot_wireframe(dw23range,wrange,var,colors='red')
ax.plot_wireframe(dw23range,wrange,-1*var,colors='red')
#ax.scatter(scatterx,scattery,scatterz)
ax.set_xlabel('dw23')
ax.set_ylabel('wness')
plt.show()


# In[125]:

raw_input()


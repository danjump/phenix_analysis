import numpy as np
#import scipy.interpolate as interpolate
from sklearn.gaussian_process import GaussianProcess
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

with open('example_a0q1.txt') as fh:
    M = np.vstack(map(float, r.split(' ')) for r in fh.read().splitlines())
r = np.linspace(-0.3, 0.3, M.shape[0])
c = np.linspace(0, 1, M.shape[1])
#r = np.linspace(0, 1, 20)
#c = np.linspace(-0.3, 0.3 100)

rr, cc = np.meshgrid(r, c)
vals = ~np.isnan(M)
uncert = np.sqrt(M[vals])
#gp = GaussianProcess(verbose=True, nugget=uncert)
gp = GaussianProcess(verbose=True,theta0=0.01, thetaL=.0001, thetaU=1., nugget=2.0)
gp.fit(X=np.column_stack([rr[vals],cc[vals]]), y=M[vals])
print 'parameters: ', gp.get_params()
print 'score: ', gp.score(np.column_stack([rr[vals],cc[vals]]), M[vals])

rr_cc_as_cols = np.column_stack([rr.flatten(), cc.flatten()])
#interpolated = gp.predict(rr_cc_as_cols).reshape(M.shape)
interpolated , predUncert = gp.predict(rr_cc_as_cols,eval_MSE=True)
interpolated = interpolated.reshape(M.shape)
predUncert = predUncert.reshape(M.shape)
upperBand = interpolated + predUncert
lowerBand = interpolated - predUncert

fig=pl.figure()
ax = fig.gca(projection='3d')
ax.plot_wireframe(rr,cc,interpolated,colors='blue')
ax.plot_wireframe(rr,cc,upperBand,colors='red')
ax.plot_wireframe(rr,cc,lowerBand,colors='red')
ax.scatter(rr,cc,M)
ax.set_xlabel('dw23')
ax.set_ylabel('wness')
pl.show()



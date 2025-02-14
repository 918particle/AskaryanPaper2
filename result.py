import numpy as np
from scipy.special import dawsn as F
from scipy.special import expi as Ei
import matplotlib.pyplot as plt

def gauss(x):
	return np.exp(-x*x)
def dDdx(x):
	return 1-2*x*F(x)
def G(x):
	return np.sqrt(np.pi)*F(x)-0.5*gauss(x)*Ei(x*x)
def plotit(X,Y,Z,title):
    plt.figure()
    cs = plt.contour(X,Y,Z,10) # contour() accepts complex values
    plt.clabel(cs, inline=1, fontsize=10) # add labels to contours
    plt.title(title)
    plt.savefig(title+'.png')
def result(x,z0,R0,E0,sigma_t):
	term1 = 4*x*x*np.abs(G(z0))*np.abs(G(z0))*np.exp(-2*x*x)
	term2 = 4*np.abs(G(z0))/np.abs(z0)*dDdx(x)*x*np.exp(-x*x)*np.cos(np.angle(z0)+np.angle(G(z0)))
	term3 = 1.0/np.abs(z0)/np.abs(z0)*dDdx(x)*dDdx(x)
	return np.sqrt(term1+term2+term3)

#x = np.arange(0.01,1,0.01)
#y = np.arange(0.01,1,0.01)
#X,Y = np.meshgrid(x,y)
#Z = np.zeros(X.shape)
#n = len(x)
#for i in range(n):
#	for j in range(n):
#		z = complex(X[i][j],Y[i][j])
#		Z[i][j] = G(z)

z0 = complex(10,0.1)
R0 = 1
E0 = 1
sigma_t = 0.5
x = np.arange(0,3,0.001)
n = len(x)
res = np.zeros(x.shape)
for i in range(n):
	res[i] = result(i,z0,R0,E0,sigma_t)

plt.plot(x,20*np.log10(res))
plt.show()
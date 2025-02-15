import numpy as np
from scipy.signal import convolve
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
def envelope(x,z0,R0,E0,sigma_t):
	x0 = z0.real
	y0 = z0.imag
	phi0 = np.angle(z0)
	abs_z_0 = np.abs(z0)
	abs_z_0_2 = abs_z_0*abs_z_0
	G0 = G(z0)
	phiG = np.angle(G0)
	abs_G_0 = np.abs(G0)
	abs_G_0_2 = abs_G_0
	phase1 = -2*x0*y0+phi0+phiG
	term1a = 16*(np.pi*np.pi)*abs_z_0_2*np.exp(-4*x*y0)*np.exp(-2*(x0*x0-y0*y0))
	term1b = 4*x*x*abs_G_0_2*np.exp(-2*x*x)
	term1c = -16*np.pi*abs_G_0*abs_z_0*np.exp(-2*x*y0)*x*np.exp(-x*x)*np.cos(2*x0*x+phase1)
	term1 = term1a+term1b+term1c
	term2 = 1/abs_z_0_2*(dDdx(x)*dDdx(x))
	term3a = dDdx(x)*(8*np.pi*np.sin(2*x*x0-2*x0*y0+2*phi0)*np.exp(-2*x*y0)*np.exp(-2*(x0*x0-y0*y0)))
	term3b = dDdx(x)*(-2*abs_G_0/abs_z_0*x*np.exp(-x*x)*np.sin(phi0+phiG))
	term3 = term3a+term3b
	return R0*E0*sigma_t*sigma_t*np.sqrt(term1+term2+term3)/np.sqrt(np.pi)/2
def signal_envelope(x,E0,sigma_t):
	return complex(-E0*np.sqrt(2)*sigma_t*x*np.exp(-x*x),-2*sigma_t/np.sqrt(2*np.pi)*dDdx(x))
def detector_envelope(t,R0,gamma,f_0):
	return R0*np.exp(-gamma*t)*complex(np.cos(2*np.pi*f_0*t),np.sin(2*np.pi*f_0*t))
def manual_envelope(t,E0,sigma_t,R0,gamma,f_0):
	return 0.5*np.abs(convolve(signal_envelope(t/np.sqrt(2)/sigma_t,E0,sigma_t),detector_envelope(t,R0,gamma,f_0),'mode=full'))

f_0 = 0.5 #GHz
gamma = 0.08 #GHz
sigma_t = 2 #ns
sigma_f = 1.0/(2*np.pi*sigma_t)
z0 = complex(f_0/(np.sqrt(2)*sigma_f),gamma/(2*np.pi*np.sqrt(2)*sigma_f))
R0 = 1
E0 = 1
t = np.arange(0.01,50,0.01)
n = len(t)
res1 = np.zeros(t.shape)
res2 = np.zeros(t.shape)
for i in range(n):
	res1[i] = envelope(t[i]/np.sqrt(2)/sigma_t,z0,R0,E0,sigma_t)
	res2[i] = manual_envelope(t[i],E0,sigma_t,R0,gamma,f_0)

#plt.plot(t,res1,'-')
plt.plot(t,res2,'-')
plt.axis([0.0, 51, -0.1, 10])
plt.show()
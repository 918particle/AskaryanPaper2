import numpy as np
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as F
from scipy.special import expi as Ei
from scipy.special import wofz as w
import matplotlib.pyplot as plt

def s(t,E0,sigma_t):
	return -E0*t*np.exp(-0.5*(t/sigma_t)*(t/sigma_t))
def wat(t,E0,R0,sigma_t,f0,gamma):
	J = complex(0,1)
	y = t/np.sqrt(2)/sigma_t
	zp = np.sqrt(2)*sigma_t*complex(gamma,2*np.pi*f0)
	zm = np.sqrt(2)*sigma_t*complex(gamma,-2*np.pi*f0)
	c = -J*(y-zp/2)
	d = -J*(y-zm/2)
	first_term = np.sqrt(np.pi)/2/np.sqrt(2)*R0*sigma_t*s(t,E0,sigma_t)*(w(c)+w(d))
	second_term = np.sqrt(2)*R0*E0*sigma_t*sigma_t*np.exp(-y*y)*(1+J*(c*w(c)+d*w(d)))
	return first_term+second_term

f0 = 0.8 #GHz
gamma = 1 #GHz
sigma_t = 0.25 #ns
R0 = 1
E0 = 1

t = np.arange(-10,10,0.01)
n = len(t)
result = np.zeros(t.shape,dtype=complex)
for i in range(n):
	result[i] = wat(t[i],E0,R0,sigma_t,f0,gamma)
plt.plot(t,np.real(result),'-',label="WAT real")
plt.plot(t,np.abs(hilbert(np.real(result))),'-',label="WAT env")
plt.legend()
plt.show()
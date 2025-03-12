import numpy as np
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as F
from scipy.special import expi as Ei
from scipy.special import wofz as w
import matplotlib.pyplot as plt

def G(z0):
	z0 = np.real(z0)
	return np.sqrt(np.pi)*F(z0)-0.5*np.exp(-z0*z0)*Ei(z0*z0)
def s(t,E0,sigma_t):
	return -E0*t*np.exp(-0.5*(t/sigma_t)*(t/sigma_t))
def r(t,R0,f0,gamma):
	return R0*np.cos(2*np.pi*f0*t)*np.exp(-gamma*t)
def ra(t,R0,f0,gamma):
	J = complex(0,1)
	return R0*np.exp(2*np.pi*J*f0*t)*np.exp(-gamma*t)
def math_conv(t,E0,R0,sigma_t,f0,gamma):
	J = complex(0,1)
	sigma_f = 1.0/(2*np.pi*sigma_t)
	p = (t/sigma_t)*(sigma_f/f0) - gamma/(2*np.pi*f0)
	q = 0.5*(sigma_f/f0)*(sigma_f/f0)
	z = complex(1.0/(2.0*np.sqrt(q)),-p/(2.0*np.sqrt(q)))
	U = np.sqrt(np.pi/(4.0*q))*np.real(w(z))
	dUdp = np.sqrt(np.pi)/(4*q)*(2/np.sqrt(np.pi) + 2*np.real(J*z*w(z)))
	first_term = np.sqrt(np.pi/2)*sigma_t*R0*s(t,E0,sigma_t)*U
	second_term = np.sqrt(np.pi/2)*E0*R0*sigma_t*sigma_t*(sigma_f/f0)*np.exp(-0.5*(t/sigma_t)*(t/sigma_t))*dUdp
	if(np.isinf(first_term) or np.isinf(second_term) or np.isnan(first_term) or np.isnan(second_term)):
		return 0.0
	else:
		return first_term+second_term
def math_env(t,E0,R0,sigma_t,f0,gamma):
	J = complex(0,1)
	sigma_f = 1.0/(2*np.pi*sigma_t)
	z0 = complex(f0/np.sqrt(2)/sigma_f,gamma/(2*np.pi*np.sqrt(2)*sigma_f))
	first_term = -E0*sigma_t*sigma_t*ra(t,R0,f0,gamma)*(1-J*np.sqrt(np.pi)*z0*w(-z0))
	second_term = E0*sigma_t*sigma_t*z0*np.exp(-z0*z0)*np.sqrt(np.pi)*ra(t,R0,f0,gamma)
	second_term += np.sqrt(2/np.pi)*G(z0)*R0*sigma_t*s(t,E0,sigma_t)
	if(t<0 or np.isinf(first_term) or np.isinf(second_term) or np.isnan(first_term) or np.isnan(second_term)):
		return 0.0
	else:
		return 0.5*np.abs(first_term+J*second_term)

f0 = 0.75 #GHz
gamma = 0.5 #GHz
sigma_t = 0.2 #ns
R0 = 1
E0 = 1
dt = 0.01
T_min = -20
T_max = 108

t = np.arange(T_min,T_max,dt)
n = len(t)
result1 = np.zeros(t.shape,dtype=complex)
result2 = np.zeros(t.shape,dtype=complex)
for i in range(n):
	result1[i] = math_conv(t[i],E0,R0,sigma_t,f0,gamma)
	result2[i] = math_env(t[i],E0,R0,sigma_t,f0,gamma)
graph1 = np.real(result1)
graph2 = np.abs(hilbert(np.real(result1)))
graph3 = np.real(result2)
graph3 *= np.max(graph2)/np.max(graph3)
plt.plot(t,graph1,'-',label="mathematical convolution")
plt.plot(t,graph2,'-',label="envelope of mathematical convolution")
plt.plot(t,graph3,'-',label="mathematical envelope")
plt.legend()
plt.show()
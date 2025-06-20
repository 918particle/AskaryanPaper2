import numpy as np
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as D
from scipy.special import wofz as w
from scipy.special import erf
import matplotlib.pyplot as plt

def dwdq(q):
	J = complex(0,1)
	return 2*J/np.sqrt(np.pi)-2*q*w(q)

def s(t,E0,sigma_t):
	return -E0*t*np.exp(-0.5*t*t/sigma_t/sigma_t)

def r(t,R0,f0,gamma):
	if(t>=0):
		return R0*np.cos(2*np.pi*f0*t)*np.exp(-2*np.pi*gamma*t)
	else:
		return 0

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
	units = R0*E0*sigma_t*sigma_t
	J = complex(0,1)
	x = t/(sigma_t*np.sqrt(2))
	z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
	q = -J*(x+z/2)
	first_term = -np.sqrt(np.pi)*units*(x*np.exp(-x*x)*w(q)+0.5*J*np.exp(-x*x)*dwdq(q))
	second_term = 0
	result = first_term+J*second_term
	if(np.isinf(result) or np.isnan(result)):
		return 0
	else:
		return 0.5*np.abs(result)

f0 = 0.33 #GHz
gamma = 0.04 #GHz
sigma_t = 1 #ns
R0 = 1
E0 = 1
dt = 0.01 #ns
T_min = -50 #ns
T_max = 50 #ns
t = np.arange(T_min,T_max,dt)
n = len(t)

J = complex(0,1)
x = t/(sigma_t*np.sqrt(2))
z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
q = -J*(x+z/2)
a=4

result1 = np.zeros(t.shape,dtype=complex)
result2 = np.zeros(t.shape,dtype=complex)
for i in range(n):
	u = t[i]
	result1[i] = math_conv(u,E0,R0,sigma_t,f0,2*np.pi*gamma)
	u = u+0.5
	u*=0.2
	result2[i] = 2*np.exp(-0.5*u*u)*(0.5*(1+erf(a*u/np.sqrt(2))))/1.6
graph1 = np.abs(hilbert(np.real(result1)))*1000
graph2 = np.real(result2)*1000
plt.plot(t,graph1,'-',label="env of math conv")
plt.plot(t,graph2,'-',label="math env")
plt.xlabel("Time [ns]")
plt.ylabel("Amplitude [mV]")
plt.xticks()
plt.yticks()
plt.legend()
plt.show()
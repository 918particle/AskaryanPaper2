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
	x = t/(sigma_t*np.sqrt(2))
	z = 2*np.pi*J*f0*np.sqrt(2)*sigma_t-2*np.pi*gamma*np.sqrt(2)*sigma_t
	k = 2*x+z
	I = np.sqrt(np.pi)/2.0*w(-J*k/2)
	dIdk = 0.5+0.25*np.sqrt(np.pi)*w(-J*k/2)
	first_term = -2*R0*E0*sigma_t*sigma_t*(x*np.exp(-x*x)*I-np.exp(-x*x)*dIdk)
	second_term = 2/np.sqrt(np.pi)*R0*E0*sigma_t*sigma_t*(F(x)-z*F(x))
	result = first_term+J*second_term
	if(np.isinf(result) or np.isnan(result)):
		return 0
	else:
		return 0.5*np.abs(result)

f0 = 0.5 #GHz
gamma = 0.9 #GHz
sigma_t = 0.4 #ns
R0 = 1
E0 = 1
dt = 0.01 #ns
T_min = -20 #ns
T_max = 120 #ns

t = np.arange(T_min,T_max,dt)
n = len(t)
result1 = np.zeros(t.shape,dtype=complex)
result2 = np.zeros(t.shape,dtype=complex)
for i in range(n):
	result1[i] = math_conv(t[i],E0,R0,sigma_t,f0,gamma)
	result2[i] = math_env(t[i],E0,R0,sigma_t,f0,gamma)
graph1 = np.real(result1)*1000
graph2 = np.abs(hilbert(np.real(result1)))*1000
graph3 = np.real(result2)*1000
#graph3 *= np.max(graph2)/np.max(graph3)
plt.plot(t,graph1,'-',label="mathematical convolution")
plt.plot(t,graph2,'-',label="envelope of mathematical convolution")
plt.plot(t,graph3,'-',label="mathematical envelope")
plt.xlabel("Time [ns]")
plt.ylabel("Amplitude [mV]")
plt.xticks()
plt.yticks()
plt.legend()
plt.show()

'''f = open("results_June5_1.dat","w")
separator = " "
for i in range(n):
	data = [t[i],graph1[i],graph2[i],graph3[i]]
	datastring = separator.join(str(x) for x in data)
	f.write(datastring+"\n")
f.close()'''
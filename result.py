import numpy as np
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as D
from scipy.special import wofz as w
import matplotlib.pyplot as plt
from scipy.integrate import quad

def shifted_laplace_transform(f,s,t0,t_max=1000):
    real_integral = quad(lambda t: np.exp(-s.real*t)*np.cos(s.imag*t)*f(t-t0),0,t_max)[0]
    imag_integral = quad(lambda t: np.exp(-s.real*t)*np.sin(s.imag*t)*f(t-t0),0,t_max)[0]
    return complex(real_integral,-imag_integral)

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
    units = R0*E0*sigma_t*sigma_t
    x = t/(sigma_t*np.sqrt(2))
    z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
    q = -J*(x+0.5*z)
    Re_wq = np.real(w(q))
    Re_dwdx = np.real(-J*dwdq(q))
    result = -np.sqrt(np.pi)*units*(x*np.exp(-x*x)*Re_wq-0.5*np.exp(-x*x)*Re_dwdx)
    if(np.isinf(result) or np.isnan(result)):
        return 0.0
    else:
        return result

def math_env(t,E0,R0,sigma_t,f0,gamma):
	units = R0*E0*sigma_t*sigma_t
	J = complex(0,1)
	x = t/(sigma_t*np.sqrt(2))
	z = (2*np.pi*J*f0-2*np.pi*gamma)*np.sqrt(2)*sigma_t
	k=-z
	q = -J*(x+z/2)
	first_term = -np.sqrt(np.pi)*units*(x*np.exp(-x*x)*w(q)+0.5*J*np.exp(-x*x)*dwdq(q))
	second_term = 2/np.sqrt(np.pi)*units*(D(x) + k*shifted_laplace_transform(D,k,x))
	if(np.isinf(first_term) or np.isinf(second_term) or np.isnan(first_term) or np.isnan(second_term)):
		return 0.0
	else:
		result = first_term+J*second_term
	return 0.5*np.abs(result)

f0 = 0.15 #GHz
gamma = 0.033 #GHz
sigma_t = 1.0 #ns
R0 = 1
E0 = 1
dt = 0.01 #ns
T_min = -100 #ns
T_max = 100 #ns
t = np.arange(T_min,T_max,dt)
n = len(t)

result1 = np.zeros(t.shape,dtype=complex)
result2 = np.zeros(t.shape,dtype=complex)
signal = np.zeros(t.shape,dtype=float)
response = np.zeros(t.shape,dtype=float)
for i in range(n):
    u = t[i]
    signal[i] = s(u,E0,sigma_t)
    response[i] = r(u,R0,f0,gamma)
    result1[i] = math_conv(u,E0,R0,sigma_t,f0,gamma)
    result2[i] = math_env(u,E0,R0,sigma_t,f0,gamma)
graph1 = np.real(result1)*1000
graph2 = np.abs(hilbert(np.real(result1)))*1000
graph3 = np.real(result2)*1000
graph4 = convolve(signal,response,"same")*dt*1000

'''plt.plot(t,graph4,'-',label="numerical convolution")
plt.plot(t,graph1,'-',label="math conv")
plt.xlabel("Time [ns]")
plt.ylabel("Amplitude [mV]")
plt.xticks()
plt.yticks()
plt.legend()
plt.show()'''

with(open('results_July7th_2.dat','w') as file):
	for i in range(n):
		write_string = str(t[i])+" "+str(graph4[i])+" "+str(graph1[i])+"\n"
		file.write(write_string)
file.close()

'''plt.plot(t,graph1,'-',label="env of math conv")
plt.plot(t,graph2,'-',label="math env")
plt.plot(t,graph3,'-',label="math conv")
plt.xlabel("Time [ns]")
plt.ylabel("Amplitude [mV]")
plt.xticks()
plt.yticks()
plt.legend()
plt.show()'''

'''with(open('results_July3rd_2.dat','w') as file):
	for i in range(n):
		write_string = str(t[i])+" "+str(graph1[i])+" "+str(graph2[i])+" "+str(graph3[i])+"\n"
		file.write(write_string)
file.close()'''
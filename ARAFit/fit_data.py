import csv
import numpy
import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as D
from scipy.special import wofz as w
from scipy.integrate import quad
from scipy.signal import correlate,correlation_lags

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
def get_csw(trace_set,index):
    ref = trace_set[index]
    n = len(ref)
    sum_trace = np.copy(ref)
    lag_array = correlation_lags(n,n,mode='full')
    for i, trace in enumerate(trace_set):
        corr = correlate(ref,trace_set[i],mode='full')/n
        best_lag = lag_array[np.argmax(corr)]
        aligned_trace = np.roll(trace,best_lag)
        sum_trace += aligned_trace
    sum_trace -= np.mean(sum_trace)
    norm = np.sqrt(np.inner(sum_trace,sum_trace))
    sum_trace/=norm
    return sum_trace

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Path to the file to open")
args = parser.parse_args()

#Data: 16 channels, 512 samples per channel, 1.5 GHz sampling rate
fs = 1.5
i=0
n = 512
data = np.zeros((7,n))
times = np.arange(0,n)/fs
#good_channels = [0,1,2,3,4,6,7] for event starting with 1915

with open(args.filename, newline='') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		data[i] = list(map(float,row))
		i+=1

csw = -1.0*get_csw(data,6)
model = np.zeros(n)
rho_max = 0
f0s = [0.150,0.155,0.160,0.165,0.170,0.175,0.180,0.185,0.190,0.195,0.200,0.205,0.210]
gammas = [0.010,0.011,0.012,0.013,0.0135,0.014,0.0145,0.015,0.016,0.017,0.018,0.019,0.020,0.021,0.022]
sigmas = [0.325,0.35,0.4,0.425,0.45,0.475,0.5,0.525,0.550]
best_sigma = 0
best_gamma = 0
best_f0 = 0
for f0 in f0s:
	for gamma in gammas:
		for sigma in sigmas:
			for i in range(n):
				model[i] = math_conv(times[i]-154,1,1,sigma,f0,gamma)
			model = model/np.sqrt(np.inner(model,model))
			rho = np.max(np.abs(correlate(model,csw)))
			if(rho>rho_max):
				best_f0=f0
				best_gamma=gamma
				best_sigma=sigma
				rho_max = rho
'''print("Best fit values")
print(best_sigma,best_f0,best_gamma,rho_max)'''
print("Data:")
for i in range(n):
	model[i] = math_conv(times[i]-156.5,1,1,best_sigma,best_f0,best_gamma)
model = model/np.sqrt(np.inner(model,model))
for i in range(n):
	print(times[i],csw[i],model[i])
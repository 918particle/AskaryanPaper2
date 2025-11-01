import csv
import numpy
import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import convolve,hilbert
from scipy.special import dawsn as D
from scipy.special import wofz as w
from scipy.integrate import quad
from scipy.signal import correlate

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

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Path to the file to open")
args = parser.parse_args()

#Data: 16 channels, 512 samples per channel, 1.5 GHz sampling rate
fs = 1.5
data = np.zeros((16,512))
i=0
times = np.arange(0,512)/fs

with open(args.filename, newline='') as csvfile:
	reader = csv.reader(csvfile, delimiter=',')
	for row in reader:
		data[i] = list(map(float,row))
		i+=1

model = np.zeros(512)
for i in range(512):
	model[i] = math_conv(times[i]-125,1,1,0.45,0.180,0.0125)
model = model/np.sqrt(np.inner(model,model))
data_n = data[3]/np.sqrt(np.inner(data[3],data[3]))
print(np.max(np.abs(correlate(model,data_n))))
plt.plot(times,model)
plt.plot(times,data_n)
plt.show()
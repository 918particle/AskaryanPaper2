import utility
import argparse
import numpy as np
from scipy.signal import correlate

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="Path to the file to open")
args = parser.parse_args()

#Data: 16 channels, 512 samples per channel, 1.5 GHz sampling rate
fs = 1.5
trace_length = 512
n_channels = 8
data = utility.load_csv_data(open(args.filename,newline=''),trace_length,n_channels)
times = np.arange(0,trace_length)/fs
good_channels = [0,3,6,7]

for k in range(8):
	if k in good_channels:
		csw = data[k]
		model = np.zeros(trace_length)
		rho_max = 0
		f0 = 0.190
		gamma = 0.014
		sigmas = np.arange(0.1,1,0.001)
		best_sigma = 0
		for sigma in sigmas:
			for i in range(trace_length):
				model[i] = utility.math_conv(times[i],1,1,sigma,f0,gamma)
			model = model/np.sqrt(np.inner(model,model))
			rho = np.max(np.abs(correlate(model,csw)))
			if(rho>rho_max):
				best_sigma=sigma
				rho_max = rho
		#print(k,best_sigma,rho_max)
	else:
		continue

csw = utility.get_csw(data,6,good_channels)

for i in range(trace_length):
	model[i] = utility.math_conv(times[i]-97,1,1,best_sigma,f0,gamma)
model = model/np.sqrt(np.inner(model,model))
for i in range(trace_length):
	print(times[i],csw[i],model[i])
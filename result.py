import numpy as np
from scipy.special import dawsn as F
from scipy.special import expi as Ei
import matplotlib.pyplot as plt

def gauss(x):
	return np.exp(-x*x)
def G(x):
	return np.sqrt(np.pi)*F(x)-0.5*gauss(x)*Ei(x*x)

x = np.arange(0,5,0.01)
plt.plot(x,G(x))
plt.plot(x,np.exp(-x)+0.185)
plt.show()

#Next: extend to complex arguments
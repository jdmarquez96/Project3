import numpy as np
import matplotlib.pyplot as plt
import math 
from scipy.optimize import curve_fit as fit

#data of the 1s upsilon and 2s upsilon
m0 = 8000 #MeV
s0 = 0.054 #MeV


Nexp = 10000
Nmeas = 1000
detector = 200




def h0(sigdetector):
        mu = np.random.normal(m0, s0) #true measurements
        x0 = np.random.normal(mu,sigdetector) #detector measurements
        return x0 



x0list = [] #generating random x cooridantes for H0
for i in range(1, Nexp):
        y0 = h0(detector)
        x0list.append(y0)
plt.figure()
plt.hist(x0list, 100, density=True) #gets x and y data list to use curve_fit
plt.title("Mass Distributions of the 1S Upsilon")
plt.xlabel("mass (GeV)")
plt.ylabel("count")
plt.savefig("massdistrib.png")
plt.close()

def gauss(x, mean, sigma):
	return (1/np.sqrt(2*np.pi*sigma**2))*np.exp(  -(x - mean)**2/(2*sigma**2))

#bin = bin0[:-1] #make sure bin and n list are same length


pullmlist = []
pullslist = []
mlist = []
slist = []
mlisterror = []
slisterror = []

for i in range(1, Nexp):
	z0list = []
	for j in range(1, Nmeas):
		y0 = h0(detector)
		z0list.append(y0)
	n, _, _ = plt.hist(z0list, 100, density = True)
	_, bin, _ = plt.hist(z0list, 100, density = True)
	binp = bin[:-1]
	param, paramcov = fit(gauss, binp, n, p0=(8000, 200))
	m = param[0]
	s = param[1]
	merr = (paramcov[0,0])
	serr = (paramcov[1,1])
	mlist.append(m)
	slist.append(s)
	mlisterror.append(merr)
	slisterror.append(serr)
	pullm = (m - m0)/merr
	pulls = (s - detector)/serr
	pullmlist.append(pullm)
	pullslist.append(pulls)
	z0list.clear()
	i += 1

plt.close()

w,_,_  = plt.hist(pullmlist, 50, density = True)
j,_,_  = plt.hist(pullslist, 50, density = True)
_, g,_ = plt.hist(pullmlist, 50, density = True)
_, k,_ = plt.hist(pullslist, 50, density = True)

param, paramcov = fit(gauss, g[:-1], w, p0 =( 0, 1))
param1, paramcov1 = fit(gauss, k[:-1], j, p0 =( 0, 1))
plt.close()

plt.figure()
plt.grid(True)
plt.hist(pullmlist, 50, density=True, label="$\\mu$ true = %.2f" %(m0))
plt.xlabel("pull values")
plt.ylabel("count")
plt.title('Pull values for $\mu$ Nexp/Nmeas 10000/100')
plt.axvline(param[0], label = '$\mu$ pull = %.2f' %(param[0]))
plt.axvline(param[0], label = '$\sigma$ pull = %.2f' %(param[1]))

plt.legend()
plt.savefig("pullmu_10000_1000_mu8000.png")
plt.show()

plt.figure()
plt.grid(True)
plt.hist(pullslist, 50, range=[-1, 1], density = True, label = "$\\sigma$ true = %.2f" %(detector))
plt.xlabel("pull values")
plt.ylabel("count")
plt.title('Pull values for $\sigma$ Nexp/Nmeas 10000/1000')
#plt.scatter(mlist,slist, yerr=slisterror, xerr=mlisterror, fmt='o')
#plt.xlabel('$\mu$ best-fit')
#plt.ylabel('$\sigma$ best-fit')
plt.axvline(0, label = '$\mu$ pull = %.2f' %(0))
plt.axvline(0, label = '$\sigma$ pull = %.2f' %(param1[1]))

plt.legend()


#plt.savefig('pullsigma_10000_1000_sigma200.png')
plt.show()


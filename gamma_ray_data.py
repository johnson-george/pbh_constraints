import math
import numpy as np
import matplotlib.pyplot as plt

###########################################################
# amalgamating the data from four places
###########################################################

f = open("igrbdata.txt", "r")
lines = f.readlines()
f.close()

def data_test(E): # temporary approximation to data, returning log10(I)	
	return -4*np.log10(E*1000)/3 - 2
logEtest = np.linspace(-6,3,100)
Etest = 10**logEtest
logItest = data_test(Etest)
#plt.plot(logEtest, logItest)

###########################################################
# HALO fits
###########################################################

def I_HALO(E): 
	Ek = E*10**6
	if Ek < 60:
		return 7.877*np.power(Ek,-0.29)*np.exp(-Ek/41.13)
	else:
		a = 0.0259*np.power(Ek/60, -5.5)
		b = 0.504*np.power(Ek/60, -1.58)
		c = 0.0288*np.power(Ek/60, -1.05)
		return a + b + c

logE = np.linspace(-5,-3,100)
logI = [np.log10(I_HALO(10**(x))) for x in logE]
#plt.plot(logE, logI)

def swifthalo(E): # the fits running down to 1 keV
	Ek = E*10**6
	C = 0.1015
	EB = 29.99
	gam1 = 1.32
	gam2 = 2.88
	den = (Ek/EB)**gam1 + (Ek/EB)**gam2
	return Ek*C/den
	
###########################################################
# FERMI LAT data
###########################################################

fermi = lines[2:27]
fermie = []
fermii = []
fermiu = []
fermitot = []
for entry in fermi:
	x = entry[:-1]
	y = [z.strip() for z in x.split(',')]
	GeVlower = float(y[0])
	GeVupper = float(y[1])
	GeV = np.exp((np.log(GeVlower) + np.log(GeVupper))/2)
	IGRBintensity = float(y[2])
	IGRBuncertainty = float(y[3])
	EGBintensity = float(y[5])
	fermie.append(GeV)
	fermii.append(IGRBintensity)
	fermiu.append(IGRBuncertainty)
	fermitot.append(EGBintensity)

#plt.plot(np.log10(fermie), np.log10(fermitot))
#plt.plot(np.log10(fermie), np.log10(fermii))

###########################################################
# EGRET data
###########################################################

egret = lines[31:42]
egrete = []
egreti = []
egretu = []
for entry in egret:
	x = entry[:-1]
	y = [z.strip() for z in x.split(',')]
	MeVlower = float(y[0])
	MeVupper = float(y[1])
	MeV = np.exp((np.log(MeVlower) + np.log(MeVupper))/2)
	deltaE = MeVupper - MeVlower
	EGRBintensity = float(y[2])
	EGRBuncertainty = float(y[3])
	egrete.append(MeV/1000)
	egreti.append(EGRBintensity*deltaE)
	egretu.append(EGRBuncertainty*deltaE)

#plt.plot(np.log10(egrete), np.log10(egreti))

###########################################################
# COMPTEL data
###########################################################

compt = lines[46:55]
compte = []
compti = []
for entry in compt:
	x = entry[:-1]
	y = [z.strip() for z in x.split(',')]
	MeV = 10**float(y[0])
	Intensity = 10**float(y[1])/MeV
	compte.append(MeV/1000)
	compti.append(Intensity)

#plt.plot(np.log10(compte), np.log10(compti))

###########################################################
# defining an interpolating function
###########################################################

egrete2 = [E for E in egrete if E < fermie[0]]
egreti2 = egreti[0:len(egrete2)]
Es = compte + egrete2 + fermie
Is = compti + egreti2 + fermii

def data(E): # actual experimental data --- returning I, not log10(I)!
	if E < Es[0]:
		return I_HALO(E)
	elif E < Es[-1]:
		j = next(i for i,f in enumerate(Es) if f>E)
		p = (np.log(E) - np.log(Es[j-1]))/(np.log(Es[j]) - np.log(Es[j-1]))
		return np.exp(np.log(Is[j-1]) + p*(np.log(Is[j]) - np.log(Is[j-1])))
	else:
		return 10**10

a = np.linspace(-6,2,1000)
energies = 10**a
b = [np.log10(data(e)) for e in energies]
c = [np.log10(e*data(e)) for e in energies]
plt.plot(a, b)
plt.show()

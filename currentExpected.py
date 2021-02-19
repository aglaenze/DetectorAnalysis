from math import *

conv = 10**3/10** (-2) 	    # conversion factor from kV/cm to V/m

gain = 2000
ampField = 26.5*conv	    # kV/cm
ibf = 0.01                  # fraction of ions that drift back
driftField = 400* 10**2     # V/cm
rate=350000                 # source rate
timeSampling=1            # time over which current is integrated, in s

qe = 1.60217662 * 10**(-19)
m = 39.948 * 1.66 * 10**(-27)   # mass or Ar ion in kg
damp = 128 * 10**(-6)
ddrift = 5 * 10**(-3)

def v(field, d):
	return 2/3*sqrt(2*qe*field/m)*sqrt(d)
 
def iExpected(field, d, N):
    return rate*timeSampling* N * qe * 2/3 * v(field, d) * d

 
iMesh=iExpected(ampField, damp, gain)
iDrift=iExpected(driftField, ddrift, ibf*gain)

print("i expected for gain " + str(gain) + " = " + str(iMesh*10**9) + " nA" )
print("i expected on drift electrode = " + str(iDrift*10**9) + " nA" )

#print("i expected for gain " + str(gain) + " = " + str(iMesh*10**9) + " nA" )


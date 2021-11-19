from math import *

# units in mm

# MGEM1
#dDrift=10.5
#dMeshTop=0.660
#dGem=0.060
#dTransfert=7.22
#dMeshBot=0.0128

# MGEM3
dDrift=13.8
dMeshTop=0.128
dGem=0.060
dTransfert=5.14
dMeshBot=0.0128


#dvGem=250   # variable à choisir
#dvMeshTop=0.15*dvGem/dGem*dMeshTop

dvMeshTop=25*400*dMeshTop/10
dvGem=dvMeshTop*dGem/dMeshTop/0.2

dvTransfert=dvGem/dGem*dTransfert/100

print('dV(mesh top) = ' + str(dvMeshTop) + ' V')
print('dV(Gem) = ' + str(dvGem) + ' V')
print('dV(transfert) = ' + str(dvTransfert) + ' V')


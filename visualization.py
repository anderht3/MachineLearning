import sys
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.special as spc
import math
from scipy.special import lpmn
import scipy.integrate as integrate
from scipy.integrate import quad
from numpy import sin, cos

# Set the number of sources and the coordinates for the input

DPI = 100
SIZE = 400

nsources = int(1297803)
nside = 8
npix = hp.nside2npix(nside)

# Coordinates and the density field f
#thetas = np.random.random(nsources) * np.pi
#phis = np.random.random(nsources) * np.pi * 2.

#fs = np.random.randn(nsources)

with open("/store/user/anderht3/cent0_default.txt") as inputFile:
    firstLine = inputFile.readline()
    lines = inputFile.readlines()
#print("read file")
events = int(firstLine)
thetas = []
phis = []

i = 0
for x in range(events):
   i += 1
   j = x+1
   #print("event ", x, " processed")
   while i < len(lines) and len(lines[i].split())< 3 and float(lines[i].split()[0]) != j:
       thetas.append(float(lines[i].split()[0]))
       phis.append(float(lines[i].split()[1]))
       i += 1
indices = hp.ang2pix(nside, thetas, phis)
hpxmap2 = np.zeros(npix, dtype = np.float)

for l in range(len(thetas)):
   hpxmap2[indices[l]] += 1.0
#print("total map made")
clFinal = np.zeros(24)
cl2Final = np.zeros(24)
clSqr = np.zeros(24)
clErr = np.zeros(24)
cl2Sqr = np.zeros(24)
cl2Err = np.zeros(24)

i = 0
for x in range(events):
   i += 1
   j = x+1
   thetas2 = []
   phis2 = []
   #print("event ", j, " processed")
   while i < len(lines) and len(lines[i].split())< 3 and float(lines[i].split()[0]) != j: 
      thetas2.append(float(lines[i].split()[0]))
      phis2.append(float(lines[i].split()[1]))
      i+=1
   
   indices2 = hp.ang2pix(nside, thetas2, phis2)
   hpxmap = np.zeros(npix, dtype = np.float)
   for k in range(len(thetas2)):
      hpxmap[indices2[k]] += 1.0
      #hpxmap[indices[k]] += npix*(1.0/len(thetas))
   hp.mollview(hpxmap, cmap = cm.jet, norm = "hist", xsize = SIZE, cbar = False, title = '')
   plt.savefig("/store/user/anderht3/ampt_visuals/def_visuals/ampt" + str(x) + ".png", dpi = DPI)
   hp.graticule()
   print(x)
   hpxmapFinal = np.zeros(npix, dtype = np.float)
   for v in range(npix):
      hpxmapFinal[v] += (hpxmap[v]/hpxmap2[v]) * (events)
   c = hp.anafast(hpxmapFinal)
   alm = hp.map2alm(hpxmapFinal, lmax = 23)
   #print(c)
   #print()
   cl2 = []
   for g in range(24):
     cl2.append(c[g] - (1.0/(2*g+1)*((abs(alm[g]))**2)))
   #print(cl2)
   for l in range(24):
     clFinal[l] = clFinal[l] + c[l]
     cl2Final[l] = cl2Final[l] + cl2[l]
     clSqr[l] = clSqr[l] + c[l]*c[l]
     cl2Sqr[l] = cl2Sqr[l] + cl2[l]*cl2[l]
for k in range (24):
   clFinal[k] = clFinal[k] / (1.0*events)
   cl2Final[k] = cl2Final[k] / (1.0*events)
   cl2Sqr[k] = cl2Sqr[k]/ (1.0*events)
   cl2Err[k] = math.sqrt((cl2Sqr[k] - cl2Final[k]*cl2Final[k])/(1.0*events))
   clSqr[k] = clSqr[k]/ (1.0*events)
   clErr[k] = math.sqrt((clSqr[k] - clFinal[k]*clFinal[k])/(1.0*events))
   print(k,clFinal[k],clErr[k])
#print()
#print()
#print(cl2Final)
print()

for k in range (24):
   print(k,cl2Final[k],cl2Err[k])

'''
with open("eventFileWCent.txt") as inputFile:
    firstLine = inputFile.readline()
    lines = inputFile.readlines()
    #print (lines[1].split()[1])

cl0 = np.zeros(24, dtype = np.float)
cl5 = np.zeros(24, dtype = np.float)
cl10 = np.zeros(24, dtype = np.float)
cl15 = np.zeros(24, dtype = np.float)
cl20 = np.zeros(24, dtype = np.float)
cl25 = np.zeros(24, dtype = np.float)
cl30 = np.zeros(24, dtype = np.float)
cl35 = np.zeros(24, dtype = np.float)

events0 = 0
events5 = 0
events10 = 0
events15 = 0
events20 = 0
events25 = 0
events30 = 0
events35 = 0

events = int(firstLine)
#print(len(lines))
i = 0
for x in range(events):
    centrality = float(lines[i].split()[1])
   
    i += 1
    j = x+1
    phis = []
    thetas = []
    while i < len(lines) and len(lines[i].split())>2 and float(lines[i].split()[0]) != j:
      #print(lines[i+1].split()[0])
      thetas.append(float(lines[i].split()[1]))
      phis.append(float(lines[i].split()[2]))
      i+=1

    indices = hp.ang2pix(nside, thetas, phis)
    hpxmap = np.zeros(npix, dtype = np.float)
    for k in range(len(thetas)):
      hpxmap[indices[k]] += 1.0
      #hpxmap[indices[k]] += npix*(1.0/len(thetas))
    hpxmapFinal = np.zeros(npix, dtype = np.float)
    for v in range(len(thetas2)):
      hpxmapFinal[indices2[v]] += (hpxmap[indices2[v]]/hpxmap2[indices2[v]]) * (npix*(1.0/len(thetas)))
    c = hp.anafast(hpxmapFinal)for k in range(len(thetas)):
      hpxmap[indices[k]] += 1.0
      #hpxmap[indices[k]] += npix*(1.0/len(thetas))
    hpxmapFinal = np.zeros(npix, dtype = np.float)
    for v in range(len(thetas2)):
      hpxmapFinal[indices2[v]] += (hpxmap[indices2[v]]/hpxmap2[indices2[v]]) * (npix*(1.0/len(thetas)))
    c = hp.anafast(hpxmapFinal)
    alm = hp.map2alm(hpxmapFinal, lmax = 23)
    cl2 = []

    alm = hp.map2alm(hpxmapFinal, lmax = 23)
    cl2 = []
    
    for g in range(len(c)):
       cl2.append(c[g] - (1.0/(2*g+1)*((abs(alm[g]))**2)))

    for z in range(len(c)):
      if (centrality < 5):
        cl0[z] = cl0[z] + cl2[z]
        events0 = events0 + 1
      if (centrality < 10 and centrality >= 5):  
        cl5[z] = cl5[z] + cl2[z]
        events5 = events5 + 1
      if (centrality < 15 and centrality >= 10):
        cl10[z] = cl10[z] + cl2[z]
        events10 = events10 + 1
      if (centrality < 20 and centrality >= 15):
        cl15[z] = cl15[z] + cl2[z]
        events15 = events15 + 1
      if (centrality < 25 and centrality >= 20):
        cl20[z] = cl20[z] + cl2[z]
        events20 = events20 + 1
      if (centrality < 30 and centrality >= 25):
        cl25[z] = cl25[z] + cl2[z]
        events25 = events25 + 1
      if (centrality < 35 and centrality >= 30):
        cl30[z] = cl30[z] + cl2[z]
        events30 = events30 + 1
      if (centrality < 40 and centrality >= 35):
        cl35[z] = cl35[z] + cl2[z]
        events35 = events35 + 1
        
for c in range(len(cl)):
    cl0[c] = cl0[c] / (1.0*events0)
    cl5[c] = cl5[c] / (1.0*events5)
    cl10[c] = cl10[c] / (1.0*events10)
    cl15[c] = cl15[c] / (1.0*events15)
    cl20[c] = cl20[c] / (1.0*events20)
    cl25[c] = cl25[c] / (1.0*events25)
    cl30[c] = cl30[c] / (1.0*events30)
    cl35[c] = cl35[c] / (1.0*events35)

#print("")
#print("")
#print("")
#print("")
#print("")
#print("")
#print("")


nsources = int(199999)
nside = 8
npix = hp.nside2npix(nside)

with open("eventFileT_iso.txt") as inputFile4:
   lines4 = inputFile4.readlines()

thetas4 = []
phis4 = []

for i in range(nsources):
   thetas4.append(float(lines4[i+1].split()[1]))
   phis4.append(float(lines4[i+1].split()[2]))

indices4 = hp.ang2pix(nside, thetas4, phis4)

hpxmap4 = np.zeros(npix, dtype = np.float)
for i in range(nsources):
    hpxmap4[indices4[i]] += npix*(1.0/nsources)

 

clm0 = np.zeros(24, dtype = np.float)
clm5 = np.zeros(24, dtype = np.float)
clm10 = np.zeros(24, dtype = np.float)
clm15 = np.zeros(24, dtype = np.float)
clm20 = np.zeros(24, dtype = np.float)
clm25 = np.zeros(24, dtype = np.float)
clm30 = np.zeros(24, dtype = np.float)
clm35 = np.zeros(24, dtype = np.float)

events0 = 0
events5 = 0
events10 = 0
events15 = 0
events20 = 0
events25 = 0
events30 = 0
events35 = 0


with open("eventFile_iso.txt") as inputFile3:
    firstLine2 = inputFile3.readline()
    lines3 = inputFile3.readlines()
    #print (lines[1].split()[1])
c2 = []
for l in range(24):
   c2.append(0)
events2 = int(firstLine2)
#print(len(lines))
i = 0
for x in range(events2):
    cent = lines3[i].split()[1]
    i += 1
    j = x+1
    phis3 = []
    thetas3 = []
    while i < len(lines3) and len(lines3[i].split()) >2  and float(lines3[i].split()[0]) != j:
      #print(lines[i+1].split()[0])
      thetas3.append(float(lines3[i].split()[1]))
      phis3.append(float(lines3[i].split()[2]))
      i+=1

    indices3 = hp.ang2pix(nside, thetas3, phis3) 
    hpxmap3 = np.zeros(npix, dtype = np.float)   
    for k in range(len(thetas3)):
      hpxmap3[indices3[k]] += npix*(1.0/len(thetas3))
    hpxmapFinal2 = np.zeros(npix, dtype = np.float)
    for v in range(len(thetas4)):
      hpxmapFinal2[indices4[v]] += hpxmap3[indices4[v]]/hpxmap4[indices4[v]]
    c2f = hp.anafast(hpxmapFinal2)



    alm2 = hp.map2alm(hpxmapFinal2, lmax = 23)
    cl2 = []

      

    for g in range(len(c2f)):
       cl2.append(c2f[g] - (1.0/(2*g+1)*((abs(alm2[g]))**2)))
    #print(alm2)
    for z in range(len(c)):
      if (centrality < 5):
        clm0[z] = clm0[z] + cl2[z]
        events0 = events0 + 1
      if (centrality < 10 and centrality >= 5):
        clm5[z] = clm5[z] + cl2[z]
        events5 = events5 + 1
      if (centrality < 15 and centrality >= 10):
        clm10[z] = clm10[z] + cl2[z]
        events10 = events10 + 1
      if (centrality < 20 and centrality >= 15):
        clm15[z] = clm15[z] + cl2[z]
        events15 = events15 + 1
      if (centrality < 25 and centrality >= 20):
        clm20[z] = clm20[z] + cl2[z]
        events20 = events20 + 1
      if (centrality < 30 and centrality >= 25):
        clm25[z] = clm25[z] + cl2[z]
        events25 = events25 + 1
      if (centrality < 35 and centrality >= 30):
        clm30[z] = clm30[z] + cl2[z]
        events30 = events30 + 1
      if (centrality < 40 and centrality >= 35):
        clm35[z] = clm35[z] + cl2[z]
        events35 = events35 + 1

for c in range(len(cl)):
    clm0[c] = clm0[c] / (1.0*events0)
    clm5[c] = clm5[c] / (1.0*events5)
    clm10[c] = clm10[c] / (1.0*events10)
    clm15[c] = clm15[c] / (1.0*events15)
    clm20[c] = clm20[c] / (1.0*events20)
    clm25[c] = clm25[c] / (1.0*events25)
    clm30[c] = clm30[c] / (1.0*events30)
    clm35[c] = clm35[c] / (1.0*events35)



for i in range(len(cl)):
   cl0[i] = cl0[i] - clm0[i]
   cl5[i] = cl5[i] - clm5[i]
   cl10[i] = cl10[i] - clm10[i]
   cl15[i] = cl15[i] - clm15[i]
   cl20[i] = cl20[i] - clm20[i]
   cl25[i] = cl25[i] - clm25[i]
   cl30[i] = cl30[i] - clm30[i]
   cl35[i] = cl35[i] - clm35[i]

'''
plt.yscale('log')
plt.plot(clFinal)
plt.savefig("powerspect_AMPT_half.png")

'''
b11 = math.sqrt((2*1+1)/(4*math.pi)*1/math.factorial(2))* -1.56727

b00 = math.sqrt((2*0+1)/(4*math.pi)*1/math.factorial(0))* 1.96962

b22 = math.sqrt((2*2+1)/(4*math.pi)*1/math.factorial(4))* 3.99862

b33 = math.sqrt((2*3+1)/(4*math.pi)*1/math.factorial(6))* -17.6705




v1 = (3/2)* cl[1]/(abs(b11)**2)*(abs(b00))**2/(4*math.pi)
v1F = math.sqrt(abs(v1))
print(v1F

v2 = (5/2)*cl0[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F0 = math.sqrt(abs(v2))
print(2.5, v2F0, 0.0)

v2 = (5/2)*cl5[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F5 = math.sqrt(abs(v2))
print(7.5, v2F5, 0.0)

v2 = (5/2)*cl10[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F10 = math.sqrt(abs(v2))
print(12.5, v2F10, 0.0)

v2 = (5/2)*cl15[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F15 = math.sqrt(abs(v2))
print(17.5, v2F15, 0.0)

v2 = (5/2)*cl20[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F20 = math.sqrt(abs(v2))
print(22.5, v2F20, 0.0)

v2 = (5/2)*cl25[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F25 = math.sqrt(abs(v2))
print(27.5, v2F25, 0.0)

v2 = (5/2)*cl30[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F30 = math.sqrt(abs(v2))
print(32.5, v2F30, 0.0)

v2 = (5/2)*cl35[2]/(abs(b22)**2)*(abs(b00))**2/(4*math.pi)
v2F35 = math.sqrt(abs(v2))
print(37.5, v2F35, 0.0)



v3 = (7/2)*cl[3]/(abs(b33)**2)*(abs(b00))**2/(4*math.pi)
v3F = math.sqrt(abs(v3))

print(v3F)



cOdd = []

for c in range(11):
   cOdd.append(cl[c*2+1])
plt.xscale('log')
plt.yscale('log')
plt.plot(cOdd)
plt.savefig("powerspect_ODD.png")
    

# Go from HEALPix coordinates to indices
indices = hp.ang2pix(nside, thetas, phis)

# Initate the map and fill it with the values
hpxmap = np.zeros(npix, dtype=np.float)
for i in range(nsources):
    #hpxmap[indices[i]] += fs[i]
    hpxmap[indices[i]] += npix*(1.0/nsources)

DPI = 100
SIZE = 400

# Inspect the map
#plt.figure(1)
'''
#map_ring = hp.pixelfunc.reorder(hpxmap, inp = 'NEST', out = 'RING')
#hp.mollview(hpxmap, xsize = SIZE)
'''
hp.mollview(hp_smoothed, cmap = cm.jet, norm = "hist", xsize = SIZE, title='Real Data smoothed')
plt.savefig("real_data_smoothed.png", dpi = DPI)
hp.graticule()
'''
'''
cl = hp.anafast(hpxmap,lmax=23)
plt.yscale('log')
axes = plt.gca()
axes.set_ylim([0.001,100])
plt.plot(cl)
plt.savefig('powerspect.png')
print(cl)
'''


'''
#plt.figure(2)
# Get the power spectrum
Cl = hp.anafast(hpxmap)
#print(Cl)
plt.plot(Cl)
plt.ylabel('C_{l}')
plt.savefig('plot_toyModel_power_spectrum.png')
'''


# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 20:38:55 2022

@author: olgso
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy_healpix import HEALPix
from astropy import units as u
#working with castles.txt
#17 sources
castles_RA=[0]*17
castles_DEC=[0]*17
f=open('castles_quasars.txt','r')
for i in range(17):
    ra_1=f.read(2)
    ra_1=float(ra_1)
    f.read(1)
    ra_2=f.read(2)
    ra_2=float(ra_2)
    f.read(1)
    ra_3=f.read(5)
    ra_3=float(ra_3)
    castles_RA[i]=15*(ra_1+ra_2/60+ra_3/3600)
    f.read(1)
    minus=f.read(1)
    dec_1=f.read(2)
    dec_1=float(dec_1)
    f.read(1)
    dec_2=f.read(2)
    dec_2=float(dec_2)
    f.read(1)
    dec_3=f.read(3)
    dec_3=float(dec_3)
    castles_DEC[i]=dec_1+dec_2/60+dec_3/3600
    if minus=='-':
        castles_DEC[i]=castles_DEC[i]*(-1)
    f.readline()
f.close()
#converting RA to -180 to 180 format
for i in range(17):
    if castles_RA[i]>180:
        castles_RA[i]=(360-castles_RA[i])*(-1)
        
hdul = fits.open('fits/String_location_ILC_direction.fits')
data = hdul[1].data
map1 = np.zeros((48, 1024),)
map2 = np.zeros((48, 1024),)
map3 = np.zeros((48, 1024),)
map = np.zeros((48 * 1024),)
long = np.zeros((48 * 1024),)
lat = np.zeros((48 * 1024),)
for i in range(48):
	map1[i, :] = np.array(data[i])[0, :]
	map2[i, :] = np.array(data[i])[1, :]
	map3[i, :] = np.array(data[i])[2, :]

hp = HEALPix(nside=64, order='ring')
for i1 in range(48):
	for i2 in range(1024):
		long_p, lat_p = hp.healpix_to_lonlat(i1 * 1024 + i2)
		if(long_p.degree <= 180.):
			long[i1 * 1024 + i2] = long_p.degree
		else:
			long[i1 * 1024 + i2] = long_p.degree - 360.
		lat[i1 * 1024 + i2] = lat_p.degree
		map[i1 * 1024 + i2] = map1[i1, i2]

in_string=[0]*17
num_in_string=0
for i in range(17):
    index=hp.lonlat_to_healpix(castles_RA[i]* u.deg,castles_DEC[i] * u.deg)
    if map[index] != 0:
        in_string[i]=1
        num_in_string=num_in_string+1
#np arrays
np_ra=np.zeros(num_in_string)
np_dec=np.zeros(num_in_string)
f=open('castles_quasars_in_strings.txt','w')
j=0;
for i in range(17):
    if in_string[i] == 1:
        np_ra[j]=castles_RA[i]
        np_dec[j]=castles_DEC[i]
        j=j+1
        f.write(str(i))
        f.write('\n')
f.close() 
 
np_castles_RA=np.array(castles_RA)
np_castles_DEC=np.array(castles_DEC) 

fig, ax = plt.subplots(3, 1)
ax[0].imshow(map1)
ax[1].imshow(map2)
ax[2].imshow(map3)
plt.show()


fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

sc = ax.scatter(long*np.pi/180., lat*np.pi/180., c = map)
#plt.colorbar(sc)

plt.subplot(projection="mollweide")
marker = '.'
marker2='+'
ms, mew = 1, 0.5
ms1,mew1 = 2,1.
ms2,mew2 =5, 1.
ms3,mew3=10,1.
plt.plot(np_castles_RA*(np.pi/180), np_castles_DEC*(np.pi/180),'ro', color='red', marker=marker, ms=ms1, mew=mew1,label='ICRF quasars')
plt.plot(np_ra*(np.pi/180), np_dec*(np.pi/180),'ro', color='white', marker=marker2, ms=ms3, mew=mew3,label='quasars in strings')
plt.show()
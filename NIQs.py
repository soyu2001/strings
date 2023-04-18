# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 21:09:20 2023

@author: olgso
"""

from astropy.io import fits
import numpy as np
from reproject import reproject_from_healpix, reproject_to_healpix
from reproject import reproject_interp
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy_healpix import HEALPix
from astropy import units as u
from photutils.segmentation import detect_sources
from astropy.wcs import WCS
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.coordinates import SkyCoord
def mollproj(ra,dec):
    th0=dec
    th=th0-(2*th0+np.sin(2*th0)-np.pi*np.sin(dec))/(2+2*np.cos(2*th0))
    while (abs(th-th0)>0.000000001):
        th0=th
        th=th0-(2*th0+np.sin(2*th0)-np.pi*np.sin(dec))/(2+2*np.cos(2*th0))
    x=2*np.sqrt(2)*(ra-180.)*np.cos(th)
    y=np.sqrt(2)*np.sin(th)
    return(x,y)
hdul = fits.open('doubles.fits')
data=hdul[1].data
cols=hdul[1].columns
RA=data['RA']
DEC=data['DEC']
names=data['Name']
target_header = fits.Header.fromstring("""
NAXIS   =                    2
NAXIS1  =                 1000
NAXIS2  =                  500
CTYPE1  = 'RA---MOL'
CRPIX1  =                500.5
CRVAL1  =                180.0
CDELT1  =                0.330
CUNIT1  = 'deg     '
CTYPE2  = 'DEC--MOL'
CRPIX2  =                250.5
CRVAL2  =                  0.0
CDELT2  =               -0.330
CUNIT2  = 'deg     '
COORDSYS= 'icrs    '
""", sep='\n')


w = WCS(target_header)

c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree)

x,y=w.world_to_pixel(c)

array, footprint = reproject_from_healpix('fits/String_location_V_direction.fits', target_header)

fig, ax = plt.subplots(figsize=(20, 10))
ax.imshow(array,interpolation='None')
array[np.isnan(array)] = 0
zerone=array
'''for i in range(500):
    for j in range(1000):
        if array[i,j]==0.:
            zerone[i,j]=0.
        else:
            zerone[i,j]=1.'''
sigma=np.std(array)
print(sigma)
            
threshold0 = array*[0]        
threshold = sigma
print(threshold)
segm = detect_sources(zerone, threshold0, npixels=10)

cmap = segm.make_cmap(seed=123)
in_string=np.zeros(len(RA))
num_in_string=0
for i in range(len(RA)):
    ix=int(x[i])
    iy=int(y[i])
    if segm.data[iy,ix] != 0:
        in_string[i]=1
        num_in_string=num_in_string+1
        #print(i,RA[i],DEC[i])
#np arrays
j=0
string_x=np.zeros(num_in_string)
string_y=np.zeros(num_in_string)
string_segm=np.zeros(num_in_string)
string_names=['nan']*num_in_string
for i in range(len(RA)):
    if in_string[i] == 1:
        string_x[j]=x[i]
        string_y[j]=y[i]
        string_segm[j]=segm.data[int(y[i]),int(x[i])]
        string_names[j]=names[i]
        j=j+1
f=open('quasars_in_strings_doubles.txt','w')
j=0;
for i in range(num_in_string):
    f.write(string_names[i])
    f.write('\n')    
f.close()
def sausage_maker(filename,len,colour):    
    def pair_file_reading(name,lines):
        ra_1=np.zeros(lines)
        ra_1_err=np.zeros(lines)
        dec_1=np.zeros(lines)
        dec_1_err=np.zeros(lines)
        ra_2=np.zeros(lines)
        ra_2_err=np.zeros(lines)
        dec_2=np.zeros(lines)
        dec_2_err=np.zeros(lines)
        yn=['nan']*lines
        f=open(name,'r')
        f.readline()
        for i in range(lines):
            s=f.readline()
            sa=s.split()
            ra_1[i]=float(sa[1])
            ra_1_err[i]=float(sa[2])
            dec_1[i]=float(sa[3])
            dec_1_err[i]=float(sa[4])
            ra_2[i]=float(sa[5])
            ra_2_err[i]=float(sa[6])
            dec_2[i]=float(sa[7])
            dec_2_err[i]=float(sa[8])
            yn[i]=sa[9]
        f.close()
        ra_1_err=ra_1_err/(1000.*3600.)
        dec_1_err=ra_1_err/(1000.*3600.)
        ra_2_err=ra_2_err/(1000.*3600.)
        dec_2_err=ra_2_err/(1000.*3600.)
        return(ra_1,ra_1_err,dec_1,dec_1_err,ra_2,ra_2_err,dec_2,dec_2_err,yn)
    def kb(ra_1,dec_1,ra_2,dec_2):
        ks=(dec_2-dec_1)/(ra_2-ra_1)
        bs=dec_1-ra_1*(dec_2-dec_1)/(ra_2-ra_1)
        k=-1./ks
        b=bs+(ra_1+ra_2)*(ks+1./ks)/2.
        return(k,b)
    ra_1,ra_1_err,dec_1,dec_1_err,ra_2,ra_2_err,dec_2,dec_2_err,yn=pair_file_reading(filename,len)
    x1=np.zeros(len)
    y1=np.zeros(len)
    x2=np.zeros(len)
    y2=np.zeros(len)
    for i in range(len):
        x1[i],y1[i]=mollproj(ra_1[i],dec_1[i])
        x2[i],y2[i]=mollproj(ra_2[i],dec_2[i])
    k=np.zeros(len)
    b=np.zeros(len)
    for i in range(len):
        k[i],b[i]=kb(ra_1[i],dec_1[i],ra_2[i],dec_2[i])
    cb = SkyCoord(ra=((ra_1+ra_2)/2.)*u.degree, dec=((dec_1+dec_2)/2.)*u.degree)
    xb,yb=w.world_to_pixel(cb)
    for i in range(len):
        base=np.linspace(xb[i],xb[i]+int(30./np.sqrt(k[i]*k[i]+1)))
        if yn[i]=='n':
            ax.scatter(base,yb[i]-k[i]*xb[i]+k[i]*base,s=0.1,c=colour)
sausage_maker('doubles_in_strings_pair_coords.txt',39,'white')
sausage_maker('NIQs_in_strings_pair_coords.txt',8,'magenta')
sausage_maker('GAIA_string_pairs_coords.txt',34,'orange')    
ci=SkyCoord('11h29m03s','+15d23m37s')
cf=SkyCoord(' 10h57m47s','+25d03m51s')
xi,yi=w.world_to_pixel(ci)
xf,yf=w.world_to_pixel(cf)
#ax.imshow(segm,cmap=cmap,interpolation='None')
ax.scatter(x,y,s=8.,c='red')
ax.scatter(xi,yi,s=30.,c='white')
ax.scatter(xf,yf,s=30.,c='white')
ax.scatter(string_x,string_y,s=8,c='white',marker='+')

'''
ra_1,ra_1_err,dec_1,dec_1_err,ra_2,ra_2_err,dec_2,dec_2_err,yn=pair_file_reading('GAIA_string_pairs_coords.txt',34)
cb1 = SkyCoord(ra=ra_1*u.degree, dec=dec_1*u.degree)
cb2 = SkyCoord(ra=ra_2*u.degree, dec=dec_2*u.degree)
x1,y1=w.world_to_pixel(cb1)
x2,y2=w.world_to_pixel(cb2)
k=np.zeros(34)
b=np.zeros(34)
for i in range(34):
    k[i],b[i]=kb(x1[i],y1[i],x2[i],y2[i])
cb = SkyCoord(ra=((ra_1+ra_2)/2.)*u.degree, dec=((dec_1+dec_2)/2.)*u.degree)
xb,yb=w.world_to_pixel(cb)
for i in range(34):
    base=np.linspace(xb[i],xb[i]+10)
    if yn[i]=='n':
        c='orange'   
        ax.scatter(base,yb[i]-k[i]*xb[i]+k[i]*base,s=0.1,c=c) 
ra_1,ra_1_err,dec_1,dec_1_err,ra_2,ra_2_err,dec_2,dec_2_err,yn=pair_file_reading('NIQs_in_strings_pair_coords.txt',8)
cb1 = SkyCoord(ra=ra_1*u.degree, dec=dec_1*u.degree)
cb2 = SkyCoord(ra=ra_2*u.degree, dec=dec_2*u.degree)
x1,y1=w.world_to_pixel(cb1)
x2,y2=w.world_to_pixel(cb2)
k=np.zeros(8)
b=np.zeros(8)
for i in range(8):
    k[i],b[i]=kb(x1[i],y1[i],x2[i],y2[i])
cb = SkyCoord(ra=((ra_1+ra_2)/2.)*u.degree, dec=((dec_1+dec_2)/2.)*u.degree)
xb,yb=w.world_to_pixel(cb)
for i in range(8):
    base=np.linspace(xb[i],xb[i]+10)
    if yn[i]=='n':
        c='magenta'   
        ax.scatter(base,yb[i]-k[i]*xb[i]+k[i]*base,s=0.1,c=c)
        '''
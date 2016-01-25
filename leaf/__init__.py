#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import numpy as np
import scipy.interpolate
import getdata
from scipy.interpolate import interp1d
import dataSpec_P5B

def trans_prosail ( N, cab, car, cbrown, cw, cm, lai, lidfa, lidfb, rsoil, psoil, \
        hspot, tts, tto, psi, typelidf):
    """A version of PROSAIL that uses transformed parameters to quasi-linearise
    the   model. See http://dx.doi.org/10.1016/j.rse.2011.12.027"""
    # Define the constants
    slai = -2.0
    skab = -100.0
    skar = -100.0
    skw =  -1./50.
    skm =  -1./100.
    # Transform the parameters to real units
    xlai = slai * np.log ( lai )
    xkab = skab * np.log ( cab )
    xkar = skar * np.log ( car )
    xkw = skw * np.log ( cw )
    xdm = skm * np.log ( cm )
    # Run the PROSAIL model
    retval = run_prosail ( N, xkab, xkar, cbrown, xkw, xdm, xlai, \
            lidfa, lidfb, rsoil, psoil, hspot, tts, tto, psi, typelidf )
    return retval

# nested stuff
'''
Define soil
'''


charsoil = '''
350 0.04,1220 0.05,1790 0.06,2300 0.07,2500 0.06'''

def soil(x,scale=None,trans=False):
    '''
    x is a dictionary
    
    ************************************
    x['spectra'] should contain spectra
    ************************************
    
    e.g.
    
    x['spectra']['dry']
    x['spectra']['char']
    
    etc.
    
    If this doesnt exist, it is loaded from the fortran data files
    and provides 'cab', 'car', 'cbrown', 'cw', 'cm'
    
    ************************************
    x['params'] should contain concentration values
    ************************************
    
     e.g.
    
    x['params']['cbrown']
    x['params']['cab']
    
    etc.
    
    '''
    # wavelength
    if not 'lamdba' in x:
        #try:
        #    x['lambda'] = eval('mod_dataspec_p5b.lambda')
        #except:
        x['lambda'] = np.arange(400.,2501.)

    x['nw'] = len(x['lambda'])

    if not 'params' in x:
        # return zero arrays of correct size if no parameters specified
        return np.zeros(x['nw']),np.zeros(x['nw']),x

    if not 'spectra' in x:
        x['spectra'] = {}

    # load defaults
    result = np.zeros(x['nw'])
    for p in x['params']:
        if not p in x['spectra']:
            if p is 'dry':
              x['spectra'][p] = mod_dataspec_p5b.rsoil1\
				-np.min(mod_dataspec_p5b.rsoil1)
              x['spectra'][p] /= x['spectra'][p].max()
            if p is 'wet':
              x['spectra'][p] = mod_dataspec_p5b.rsoil2\
				-np.min(mod_dataspec_p5b.rsoil2)
              x['spectra'][p] /= x['spectra'][p].max()
            if p is 'char':
              cchar = np.array([np.array(i.split()).astype(float) for i in charsoil.split(',')]).T
              x['spectra'][p] = \
                scipy.interpolate.interp1d(cchar[0],cchar[1])(x['lambda'])
        try:
                result += x['params'][p] * x['spectra'][p]
        except:
                pass
    return result,x



import numpy as np
from getdata import getdata
import os
import sys

'''
Define leaf absorbing constituents
'''

class Leaf():
  def __init__(self,theta=40.0,store=False,verbose=False,TINY=1e-20):
    self.verbose=verbose
    self.theta = theta
    self.store=store
    self.TINY = TINY
    self.errors = []
    self.internal_db = os.path.abspath(os.path.realpath(__file__))\
		.replace('__init__.pyc','').replace('__init__.py','')
    if verbose: print self.internal_db

    self.db = getdata(self.internal_db)
    self.data = self.db.copy()
    self.nw = self.db['refractive'].shape

  def error(self,msg):
    print msg
    sys.exit(0)

  def verbose_level(self,level=0):
    self.verbose = level

  def add_spectra(self,s,conc=None):
    '''
    s is a dictionary 

    It can contain:

    'N' : number of leaf layers
    'lambda' : new spectral sampling (e.g. vis/nir only)
               You can change the spectral sampling
               by putting a new lambda in s
               then all existing data transformed to this
    'theta'  : new angle for transmission calculation

    also, arrays:
    'xx' : Either: array of shape (2,nl_) specifying
                   wavelength (nm) and absorption coefficients for 
                   term 'k_xx'. The data are interpolated
                   from the specified wavelengths (column 0)
                   to those in self.data['lambda'] 
           Or:     array of shape(nl,) of absorption
                   coefficients as a function of self.data['lambda']

    Note: specifying s['lambda'] is applied before loading any 'xx' terms.
  
    Note: interp1d used, which may have bounds fail, so should be made safe.

    '''
    if 'N' in s:
      self.data['N'] = s['N']

    if 'lambda' in s:
      for k in self.data:
        if not (k == 'lambda' or k == 'N' or k == 'theta'):
          f = interp1d(self.data['lambda'],self.data[k])
          self.data[k] = f(s['lambda'])
      self.data['lambda'] = s['lambda'] 
    
    for k in s:
      if self.verbose: print 'considering',k
      if len(k) > 2 and (k[0] == 'k' or k[0] == 'K'):
        k0 = k
        k = k[2:]
        kk = 'k_'+k
        if self.verbose:
          print 'trying',k,'for',kk
        try:
          f = interp1d(s[k0][0],s[k0][1])
          self.data[kk] = f(s['lambda'])
        except:
          if s[k].shape == self.data['lambda'].shape:
            self.data[kk] = s[k0]
          else:
            err = 'error in spectra specification in Leaf.add_spectra() for '+k+' '+kk
            self.errors.append(err)
            if self.verbose: print err
    # update
    self.nw = self.data['refractive'].shape

  def getk(self,conc):
    '''
    Interpret concentration data dictionary
    into k
    '''
    if 'N' in conc:
      self.N = conc['N']
    else:
      self.N = 1.0

    self.conc = {}
    kvalue = np.zeros(self.nw) 
    for k in conc.keys():
      # find spectra in database
      kk = 'k_' + k
      if kk in self.data:
        if self.verbose:
          print 'found',kk,'from',k,'at',conc[k]
        self.conc[kk] = conc[k]
        kvalue += conc[k] * self.data[kk]
    # it will cause issues if 0, so let it be tiny
    kvalue[np.where(kvalue<self.TINY)] = self.TINY
    return kvalue/self.N

  def rt(self,conc):
    '''
    calculate reflectance and transmittance
    for given concentrations

    '''
    self.t1 = self.tav_abs(90.)
    if 'theta' in conc:
      self.t2 = self.tav_abs(conc['theta'])
    else:
      self.t2 = self.tav_abs(self.theta)
    # should log which was used

    self.conc = conc 
    self.k   = self.getk(conc)
    # note self.N might not be set until
    # self.getk() is called
    N = self.N

    self.tau = np.zeros_like(self.k)
    # upper limit
    ww = np.where(self.k >= 85)[0]
    if len(ww): 
      self.tau[ww] = 0.
    # lower limit
    ww = np.where(self.k <= 4)[0]
    if len(ww):
      xx=0.5*self.k[ww]-1.0
      yy=(((((((((((((((-3.60311230482612224e-13 \
            *xx+3.46348526554087424e-12)*xx-2.99627399604128973e-11) \
            *xx+2.57747807106988589e-10)*xx-2.09330568435488303e-9) \
            *xx+1.59501329936987818e-8)*xx-1.13717900285428895e-7) \
            *xx+7.55292885309152956e-7)*xx-4.64980751480619431e-6) \
            *xx+2.63830365675408129e-5)*xx-1.37089870978830576e-4) \
            *xx+6.47686503728103400e-4)*xx-2.76060141343627983e-3) \
            *xx+1.05306034687449505e-2)*xx-3.57191348753631956e-2) \
            *xx+1.07774527938978692e-1)*xx-2.96997075145080963e-1
      yy=(yy*xx+8.64664716763387311e-1)*xx+7.42047691268006429e-1
      yy=yy-np.log(self.k[ww])
      self.tau[ww] = (1.0-self.k[ww])*np.exp(-self.k[ww])+self.k[ww]**2*yy

    ww = np.where((self.k > 4) * (self.k <= 85))[0]
    if len(ww):
      xx=14.5/(self.k[ww]+3.25)-1.0
      yy=(((((((((((((((-1.62806570868460749e-12 \
                *xx-8.95400579318284288e-13)*xx-4.08352702838151578e-12) \
                *xx-1.45132988248537498e-11)*xx-8.35086918940757852e-11) \
                *xx-2.13638678953766289e-10)*xx-1.10302431467069770e-9) \
                *xx-3.67128915633455484e-9)*xx-1.66980544304104726e-8) \
                *xx-6.11774386401295125e-8)*xx-2.70306163610271497e-7) \
                *xx-1.05565006992891261e-6)*xx-4.72090467203711484e-6) \
                *xx-1.95076375089955937e-5)*xx-9.16450482931221453e-5) \
                *xx-4.05892130452128677e-4)*xx-2.14213055000334718e-3
      yy=((yy*xx-1.06374875116569657e-2)*xx-8.50699154984571871e-2)*xx+\
                9.23755307807784058e-1
      yy=np.exp(-self.k[ww])*yy/self.k[ww]
      self.tau[ww]=(1.0-self.k[ww])*np.exp(-self.k[ww])+self.k[ww]**2*yy

    tau = self.tau
    refr2 = self.data['refractive']**2
    x1=1-self.t1
    x2=self.t1**2*tau**2*(refr2-self.t1)
    x3=self.t1**2*tau*refr2
    x4=refr2*refr2-tau**2*(refr2-self.t1)**2
    x5=self.t2/self.t1
    x6=x5*(self.t1-1)+1-self.t2
    r=x1+x2/x4
    t=x3/x4
    ra=x5*r+x6
    ta=x5*t

    # to store if needed
    if self.store:
      self.c1=x6
      self.ra=x5*ra+x6
      self.ta=x5*ta
      self.c2=x1
      self.r1=x1+x2/x4
      self.t1=x3/x4
 
    '''
    reflectance and transmittance of N layers
    Stokes G.G. (1862), On the intensity of the light reflected from or transmitted
    through a pile of plates, Proceedings of the Royal Society of London, 11:545-556.
    ''' 
    r[r<self.TINY]=self.TINY
    delta=(t**2-r**2-1)**2-4*r**2
    delta[delta<0]=0.
    beta=(1+r**2-t**2-np.sqrt(delta))/(2*r)
    va=(1+r**2-t**2+np.sqrt(delta))/(2*r)
    va[va<self.TINY]=self.TINY
    beta[np.abs(beta-r)<self.TINY]+=self.TINY
    den = (va*(beta-r))
    den[den<self.TINY]=self.TINY
    ss = beta*(va-r)/den
    ss[ss<0.] = 0.
    vb=np.sqrt(ss)
    s1=ra*(va*vb**(N-1)-va**(-1)*vb**(-(N-1)))+(ta*t-ra*r)*(vb**(N-1)-vb**(-(N-1)))
    s2=ta*(va-va**(-1))
    s3=va*vb**(N-1)-va**(-1)*vb**(-(N-1))-r*(vb**(N-1)-vb**(-(N-1)))
    self.r=s1/s3
    self.t=s2/s3
 
  '''
  tau average
  '''
  def tav_abs(self,theta):
    '''
    average transmittance for given refractive index

    computation of the average transmittivity at the leaf surface within a given
    solid angle. teta is the incidence solid angle (in radian). The average angle
    that works in most cases is 40deg*pi/180. ref is the refaction index.
    ********************************************************************************
    Stern F. (1964), Transmission of isotropic radiation across an interface between
    two dielectrics, Applied Optics, 3:111-113.
    Allen W.A. (1973), Transmission of isotropic light across a dielectric surface in
    two and three dimensions, Journal of the Optical Society of America, 63:664-666.
    ********************************************************************************^
    version 5.02 (25 July 2011)


    '''
    refr = self.data['refractive']
    thetarad=np.pi*theta/180.
    if theta == 0:
        res=4.*refr/(refr+1.)**2
        return res

    refr2=refr*refr
    ax=(refr+1.)**2/2.
    bx=-(refr2-1.)**2/4.

    if thetarad == np.pi/2.:
        b1=0.
    else:
        b1=np.sqrt((np.sin(thetarad)**2-(refr2+1.)/2.)**2+bx)

    b2=np.sin(thetarad)**2-(refr2+1.)/2.
    b0=b1-b2
    ts=(bx**2/(6.*b0**3)+bx/b0-b0/2.)-(bx**2/(6.*ax**3)+bx/ax-ax/2.)
    tp1=-2.*refr2*(b0-ax)/(refr2+1.)**2
    tp2=-2.*refr2*(refr2+1.)*np.log(b0/ax)/(refr2-1.)**2
    tp3=refr2*(1./b0-1./ax)/2.
    tp4=16.*refr2**2*(refr2**2+1.)*np.log((2.*(refr2+1.)*b0-(refr2-1.)**2)/ \
            (2.*(refr2+1.)*ax-(refr2-1.)**2))/((refr2+1.)**3*(refr2-1.)**2)
    tp5=16.*refr2**3*(1./(2.*(refr2+1.)*b0-((refr2-1.)**2))-1./(2.*(refr2+1.) \
            *ax-(refr2-1.)**2))/(refr2+1.)**3
    tp=tp1+tp2+tp3+tp4+tp5
    res=(ts+tp)/(2.*np.sin(thetarad)**2)

    return res

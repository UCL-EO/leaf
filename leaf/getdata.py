#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import numpy as np
from os import sep

def getdata(here,dfile='dataSpec_P5B.py'):
  '''
  grab spectral datasets from prospect file
  '''
  try:
    lines = open(here + sep + dfile).readlines()
  except:
    try:
      lines = dataSpec_P5B.split('\n')
    except:
      self.error('unable to find database information in %s'%dfile)
      self.exit()
  data = {}
  for i in xrange(len(lines)):
    if lines[i].find('nw=')>0:
        nw = int(lines[i].split('nw=')[1])
    if lines[i].find('DATA')>0:
        # get the name, eg lambda
        term = lines[i].split('(')[1]
        # get the indices
        n = np.array(lines[i].split('i=')[1].split(')')[0].split(',')).astype(int) - 1
        # first time
        if n[0] == 0:
          data[term] = np.zeros(nw)
        # get next lines until find /
        liner = lines[i].split('i=')[1].split(')')[1].strip()
        while liner[-1] == '&':
            i += 1
            liner += lines[i].strip()
        try:
            liner = np.array(liner.replace('&','')[1:-1].split(',')).astype(float)
        except:
            liner = 0.
        # now load
        data[term][n[0]:n[1]+1] = liner
  return data

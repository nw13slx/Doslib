#!/usr/bin/env python

class atom_data(object):
  natom=None
  ntype=None
  symbol=None
  species=None
  positions=None
  constraints=None
  boundary=None
  ngh_list=None
  ngh_id=None

class dos_data(object):
  nedos = None
  efermi = None

  Xenergy = None
  dos0 = None

  spin=True
  nsites=None
  ncols=None
  dos_par=None
  loc_down = None 
  loc_up = None
  fermiN = None
  par_element=[-1,-1,-1,-1]
  par_orbital=None
  perspecies=None
  d_t2g=None
  d_eg=None
  tot=None
  partial=None


#!/usr/bin/env python

#usage: split_dos_atom argument-name argument
#argument name, parameter
# partial orbital: "s, p, d, f", # of element; i.e. p 1 f 0
# input element species: atom.ntype, [# of types, type 0 #, type 1 #, ...]
#  this setting overrides the one in poscar/contcar
# write_dos0 yes/no
# write_pdos yes/no
# center_ef yes/no
# zoomin emin emax

class control_knob:

  par_symbol=['s','p','d','f']

  def __init__(self,argv,atom,dos):
    if not argv:
      print "no input argument,use default value"

    self.name=""
    self.path=""
    self.anados=True

    self.doscar=""
    self.poscar=""

    self.run_pdos=True
    self.write_dos0=True
    self.write_pdos=False
    self.peratom=False
    self.perspecies=False
    self.whole_range=False
    self.center_ef=False

    self.zoom_in=False
    self.zoom_emax=None
    self.zoom_emin=None
    self.energyshift=0.0
    self.end=None
    self.atom=atom
    self.dos=dos

    if argv:
      for i, u in argv.__dict__.items():
        if (i=='zoomin'):
          if (u!=None):
            self.zoom_in=True
            self.zoom_emin=argv.zoomin[0]
            self.zoom_emax=argv.zoomin[1]
        elif (i in self.par_symbol):
          j=self.par_symbol.index(i)
          self.dos.par_element[j]=u
        elif (( i =="atom_ntype") and (u!=None)):
          self.atom.ntype=len(u)
          self.atom.species = []
          for j in range(atom.ntype):
            self.atom.species += ([j]*int(u[j]))
          self.atom.natom=len(self.atom.species)
        else:
          for j, v in self.__dict__.items():
            if ((i == j) and (type(u)==type(v))):  #I changed it
              #print "read parameter",i,"value", argv.__dict__[i]
              self.__dict__[j] = argv.__dict__[i]

      if  ((self.perspecies==False) and all(x==-1 for x in self.dos.par_element) and (self.peratom==False) and (self.write_pdos==False)):
        self.run_pdos=False
      else:
        self.run_pdos=True
      if (self.run_pdos==False):
        self.write_pdos=False
      if (self.write_pdos==True):
        self.perspecies=True


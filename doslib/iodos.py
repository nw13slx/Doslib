#data storage

#!/usr/bin/env python
import numpy as np
import os

class iodos:
  '''this is the class that take care of the input and output of doscar '''
  #virtual function
  tally=None

  par_symbol=['s','p','d','f']
  def _auto_terminate(self):
    if self.end:
        raise SystemExit('END: %s' % self.end)

  def __init__(self,control,plot):
    self.start=control.anados
    self.end=control.end
    self._auto_terminate()
    self._cont=control
    self.atom=control.atom
    self.dos=control.dos
    self.plot=plot

    self.end=None

    self._posf=None
    self._dosf=None

    self._find_file()

  def _find_file(self):
    '''find out atomic coordinations and DOSCAR; the priority is command line input "atom.ntype" > "poscar" input>CONTCAR >POSCAR , and "doscar" input >DOSCAR) '''
    control=self._cont
    atom=self.atom
    dos=self.dos
    if atom.ntype:
      atom.natom=len(atom.species)
    elif (os.path.isfile(control.path+control.poscar) and control.poscar):
      self._posf=open(control.path+control.poscar)
    elif os.path.isfile(control.path+"CONTCAR"):
      self._posf =open(control.path+ "CONTCAR")
    elif os.path.isfile(control.path+"POSCAR") :
      self._posf =open(control.path+"POSCAR")
    else:
      self.end="cannot find poscar@"+control.path

    if (self.start):
      if os.path.isfile(control.path+control.doscar) and control.doscar:
        self._dosf=open(control.path+control.doscar)
      elif os.path.isfile(control.path+"DOSCAR"):
        self._dosf =open(control.path+"DOSCAR")
      else:
        self.end="cannot find doscar@"+control.path
      self._auto_terminate()

  def read_poscar(self):
    '''read POSCAR or CONTCAR, with or without constraints '''
    if (self._posf):
      atom=self.atom
      lines=self._posf.readlines()
      scale=float(lines[1])
      boundary=[map(float,lines[2].strip().split())]
      boundary+=[map(float,lines[3].strip().split())]
      boundary+=[map(float,lines[4].strip().split())]
      atom.boundary=np.array(boundary)*scale
      initial=lines[5].strip().split()[0][0].lower()
      if atom.species is None:
        if ((initial>='a') and (initial<='z')):
          symbol=lines[5].strip().split()
          number=map(int,lines[6].strip().split())
          ntype=len(number)
          line_n=7
        else:
          number=map(int,lines[5].strip().split())
          ntype=len(number)
          symbol=range(ntype)
          line_n=6
      if atom.species is None:
        atom.species=[]
        for i in range(len(number)):
          atom.species += ([i]*number[i])
        atom.natom=len(atom.species)
        atom.symbol=symbol
        atom.ntype=ntype
      #if there's a line of "selective dynamics'
      #read the constraints
      initial=lines[line_n].strip().split()[0][0]
      if ('s'==initial.lower()):
          starting=line_n+2
          a=np.array([line.strip().split() for line in lines[starting:(starting+atom.natom)]])
          atom.positions=a[:,0:3].astype(np.float)
          atom.constraints=a[:,3:6]
      else:
          starting=line_n+1
          atom.positions = np.loadtxt(lines[starting:(starting+atom.natom)])
      if ('d'==(lines[starting-1].strip().split()[0][0]).lower()):
          atom.positions = atom.positions.dot(boundary)
      self._posf.close()

  def read_tot_dosfile(self):
    if (not self.start):
      return
    '''read the overall dosfile'''
    control=self._cont
    atom=self.atom
    dos=self.dos

    if not self._dosf:
        self.end="doscar is not open"
    self._auto_terminate()

    line = self._dosf.readline()
    doscaratom=int(line.strip().split()[0])
    if (doscaratom != atom.natom):
        print "WARNING: doscar atom number is different from poscar %d %d"%(doscaratom,atom.natom)
        if (doscaratom <atom.natom):
          atom.natom=doscaratom
    self._auto_terminate()

    #read the header
    for i in range(5):
      line = self._dosf.readline()
    split=line.strip().split()
    dos.nedos = int(split[2])
    dos.efermi = float(split[3])

    dos.Xenergy = np.zeros(dos.nedos)
    #read the DOS0 file
    chunck = []
    for n in xrange(dos.nedos):
        chunck.append(self._dosf.readline())
    data = np.loadtxt(chunck)
    dos.Xenergy = data[:,0]
    if (len(data[0,:])==5):
      dos.spin=True
      dos.dos0 = data[:,1:3]
      dos.dos0[:,1]=-dos.dos0[:,1]
    else:
      dos.spin=False
      dos.dos0 = data[:,1]

    if (control.center_ef==True):
      dos.Xenergy -= dos.efermi
      dos.loc_down = np.argmin(abs(dos.Xenergy+7))
      dos.loc_up = np.argmin(abs(dos.Xenergy-7))
      dos.fermiN = np.argmin(abs(dos.Xenergy))
    elif (control.energyshift!=0):
      print "shift energy by ", control.energyshift
      dos.fermiN = np.argmin(abs(dos.Xenergy-dos.efermi))
      dos.Xenergy -= control.energyshift
    else:
      dos.fermiN = np.argmin(abs(dos.Xenergy-dos.efermi))
    if (control.zoom_in==True):
      dos.loc_down = np.argmin(abs(dos.Xenergy-control.zoom_emin))
      dos.loc_up = np.argmin(abs(dos.Xenergy-control.zoom_emax))
      print "dos up down limit",dos.loc_down,dos.loc_up
    elif (control.whole_range==True):
      dos.loc_down = 0
      dos.loc_up = dos.nedos

  def write_tot_dosfile(self):
    if (not self.start):
      return
    control=self._cont
    dos=self.dos
    matrix = np.hstack([dos.Xenergy.reshape([dos.nedos,1]),dos.dos0])
    np.savetxt(control.name+'DOSall.gz',matrix,fmt='%15.8f')
    del matrix

  def write_pdosfile(self):
    if (not self.start):
      return
    if not self._cont.run_pdos:
       return 0
    control=self._cont
    dos=self.dos
    if dos.spin:
      matrix = np.hstack([dos.Xenergy.reshape([dos.nedos,1]),dos.perspecies[:,:,0],dos.perspecies[:,:,1]])
    else:
      matrix = np.hstack([dos.Xenergy.reshape([dos.nedos,1]),dos.perspecies[:,:,0]])
    np.savetxt(control.name+'DOS-perspecies.gz',matrix,fmt='%15.8f')
    del matrix
    matrix = dos.Xenergy.reshape([dos.nedos,1])
    for orbital in range(4):
      if (dos.par_orbital[orbital] is not None):
          matrix = np.hstack([matrix,dos.par_orbital[orbital]])
    np.savetxt(control.name+'DOS-parorbital.gz',matrix,fmt='%15.8f')
    del matrix

  def write_peratom(self,atomi):
    if (not self.start):
      return
    if not self._cont.run_pdos:
       return 0
    control=self._cont
    dos=self.dos
    tot=dos.tot
    matrix = np.hstack([dos.Xenergy.reshape([dos.nedos,1]),tot])
    np.savetxt(control.name+'DOS-atom'+str(atomi)+'.gz',matrix,fmt='%15.8f')
    del matrix

  def delete_tot_dosfile(self):
    self.dos.dos0=None

  def delete_pdosfile(self):
    if not self._cont.run_pdos:
       return 0
    self.dos.par_orbital=None
    self.dos.perspecies=None

  def read_pdos(self):
    if (not self.start):
      return
    if not self._cont.run_pdos:
       return 0

    dosf=self._dosf
    dos=self.dos
    atom=self.atom
    par_element=dos.par_element
    peratom=self._cont.peratom
    perspecies=self._cont.perspecies

    start = dosf.tell()
    line = dosf.readline()
    if not line:
      self.end="DOSCAR does not contain pdos information"
    self._auto_terminate()

    line = dosf.readline().strip().split()
    ncols=int(len(line))
    species=atom.species

    dosf.seek(start)
    if (ncols==3 or ncols==9 or ncols==19 or ncols==33):
       dos.spin=True
       dos.ncols=ncols
       dos.nsites = (ncols -1)/2
       self.tally=self.tally_spin
       nspin=2
    elif (ncols==2 or ncols==5 or ncols==10 or ncols==17):
       dos.spin=False
       dos.ncols=ncols
       dos.nsites=ncols-1
       self.tally=self.tally_nospin
       nspin=1
    else:
        seld.end="weird DOSCAR, I cannot read it..."
    self._auto_terminate()


    for i in range(4):
      if (par_element[i]>=0) and (dos.nsites<(2*i+1)):
        self.end="nsites is smaller than the required orbital"
      if par_element[i]>=atom.ntype:
        self.end="the type # for orbital %s is higher then the total types of element"%par_symbol[i]
    self._auto_terminate()

    if (par_element[2]>=0):
      dos.d_t2g=np.zeros((dos.nedos,nspin))
      dos.d_eg=np.zeros((dos.nedos,nspin))
    else:
      dos.d_t2g=None
      dos.d_eg=None
    dos.par_orbital=[] #store the partial s, p, d, f orbital
    for i in range(4):
      if par_element[i]>=0:
        dos.par_orbital+=[np.zeros((dos.nedos,nspin))]
      else:
        dos.par_orbital+=[None]

    if perspecies:
      if dos.spin:
        dos.perspecies=np.zeros((dos.nedos,atom.ntype,nspin))
      else:
        dos.perspecies=np.zeros((dos.nedos,atom.ntype,1))
    if peratom:
      self.plot.atom_start()
      self.atom_start()
    for atomi in xrange(atom.natom):
      #skip the first line, and loop over dos.nedos
      tot,partial=self.read_atomDOS()
      dos.tot=tot
      dos.partial=partial
      #tally the dos per orbital and per species
      for orbital in range(4):
        if (species[atomi] == par_element[orbital]):
          dos.par_orbital[orbital]+=partial[:,nspin*orbital:nspin*orbital+nspin]
      if (species[atomi] == par_element[2]):
        dos.d_t2g+=partial[:,nspin*4:nspin*4+nspin]
        dos.d_eg+=partial[:,nspin*5:nspin*5+nspin]
      if perspecies:
        if dos.spin:
          print "atomi",atomi,"species",species[atomi],"element",np.min(tot[:,0]),np.max(tot[:,0])
          dos.perspecies[:,species[atomi],0] += tot[:,0]
          dos.perspecies[:,species[atomi],1] += tot[:,1]
        else:
          dos.perspecies[:,species[atomi],0] += tot[:,0]
      if peratom:
        self.plot.atom_during(atomi)
        self.atom_during(atomi)
      del tot,partial
      #print data and png per atom
#      #if (peratom==True):
#      #  matrix = np.hstack([xEnergy.reshape([len(xEnergy),1]),partial])
#      #  header_line="#"+name+"atom"+str(atomi)+"\n#E su sd pu pd du dd fu fd dt2gu dt2gd degu degd"
#      #  np.savetxt(name+'DOS'+str(atomi)+".dat",matrix,fmt='%15.8f',header=header_line)

    if peratom:
      self.plot.atom_end()
      self.atom_end()

  def atom_start(self):
      pass
  def atom_during(self,atomi):
    control=self._cont
    if control.write_pdos :
      self.write_peratom(atomi)

  def atom_end(self):
      pass

  def read_atomDOS(self):
      if (not self.start):
        return
      dos=self.dos
      atom=self.atom
      '''read in the dos data for each atom. I put it as a separate function, in case I have some separated data in the future '''
      pdosf=self._dosf
      line = pdosf.readline()
      chunck=[]
      for n in xrange(dos.nedos):
        chunck.append(pdosf.readline())
      element = np.loadtxt(chunck)
      return self.tally(element)

  def tally_nospin(self,element):
      dos=self.dos
      atom=self.atom
      nsites = dos.nsites
      ncols=dos.ncols
      tot=np.zeros((dos.nedos,1))
      partial=np.zeros((dos.nedos,6)) #s,p,d,f,dt2g,deg
      for site in range(nsites):
        tot[:,0] += element[:,site+1]
      #collection, s, p, d, f, dt2g, deg
      partial[:,0]=element[:,1] #s
      if (ncols>4): #p
        partial[:,1]=element[:,2] + element[:,3] + element[:,4]
      if (ncols>9): #d
        partial[:,2]=element[:,5] + element[:,6] + element[:,7] + element[:,8] + element[:,9]
        partial[:,3]=element[:,5] + element[:,6] + element[:,8] #dt2g
        partial[:,4]=element[:,7] + element[:,9] #deg
      if (ncols>16): #f
        partial[:,5]=element[:,10] + element[:,11] + element[:,12] + element[:,13] + element[:,14]+element[:,15]+element[:,16]
      return tot,partial

  def tally_spin(self,element):
      dos=self.dos
      atom=self.atom
      nsites = dos.nsites
      ncols=dos.ncols
      tot=np.zeros((dos.nedos,2))
      partial=np.zeros((dos.nedos,12)) #s,p,d,f,dt2g,deg
      for site in range(nsites):
        tot[:,0] += element[:,site*2+1]
        tot[:,1] -= element[:,site*2+2]
      #collection, s, p, d, f, dt2g, deg
      partial[:,0]=element[:,1] #s
      partial[:,1]-=element[:,2] #s
      if (ncols>8): #p
        partial[:,2]=element[:,3] + element[:,5] + element[:,7]
        partial[:,3]-=element[:,4] + element[:,6] + element[:,8]
      if (ncols>18): #d
        partial[:,4]=element[:,9] + element[:,11] + element[:,13] + element[:,15] + element[:,17]
        partial[:,5]-=element[:,10] + element[:,12] + element[:,14] + element[:,16] + element[:,18]
        partial[:,8]=element[:,9] + element[:,11] + element[:,15] #dt2g
        partial[:,9]-=element[:,10] + element[:,12] + element[:,16] #dt2g
        partial[:,10]=element[:,13] + element[:,17] #deg
        partial[:,11]-=element[:,14] + element[:,18] #deg
      if (ncols>32): #f
        partial[:,6]=element[:,19] + element[:,21] + element[:,23] + element[:,25] + element[:,27]+element[:,29]+element[:,31]
        partial[:,7]-=(element[:,20] + element[:,22] + element[:,24] + element[:,26] + element[:,28]+element[:,30]+element[:,32])
      return tot,partial

def dumpclean(obj):
  '''this source code comes from http://stackoverflow.com/questions/15785719/how-to-print-a-dictionary-line-by-line-in-python'''
  if type(obj) == dict:
    for k, v in obj.items():
      if type(v) == list:  #I changed it
          print k," : ",v  #I changed it
      elif hasattr(v, '__iter__'):
         print k
         dumpclean(v)
      else:
         if v:
           print '%s : %s' % (k, v)
  else:
    print obj



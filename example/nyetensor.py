import os, sys, inspect
import operator
 # realpath() will make your script run, even if you symlink it :)
#cmd_folder = "/home/leleslx/basil/VASP/dos-lib"
#if cmd_folder not in sys.path:
#  sys.path.insert(0, cmd_folder)

from doslib.control import *
from doslib.iodos import *
from doslib.data import *
from doslib.plot import *
from doslib.analysis import *
from doslib.control import *
from doslib.iodos import *
from doslib.data import *
from doslib.plot import *
from doslib.analysis import *

import sys


NN=4
NNN=[2.6,2.85]
a=1/np.sqrt(2.)
b=0.5
latt=5.49
Ceref=np.array([[a,0,b],[0,-a,b],[-a,0,b],[0,a,b],[a,0,-b],[0,-a,-b],[-a,0,-b],[0,a,-b]])*latt
Ceref_n=Ceref
a=np.sqrt(2.)/4.
O1ref=np.array([[0,a,b],[0,-a,b],[a,0,-b],[-a,0,-b]])*latt
O1ref_n=O1ref
O2ref=np.array([[a,0,b],[-a,0,b],[0,a,-b],[0,-a,-b]])*latt
O2ref_n=O2ref
Ceref_r=np.zeros(len(Ceref))
O1ref_r=np.zeros(len(O1ref))
O2ref_r=np.zeros(len(O2ref))
DIST=0
DX=1
ANGLE=2
ID=3

for i in range(len(Ceref)):
  Ceref_r[i]=np.linalg.norm(Ceref[i,:])
  Ceref_n[i,:]/=Ceref_r[i]
for i in range(len(O1ref)):
  O1ref_r[i]=np.linalg.norm(O1ref[i,:])
  O1ref_n[i,:]/=O1ref_r[i]
  O2ref_r[i]=np.linalg.norm(O2ref[i,:])
  O2ref_n[i,:]/=O2ref_r[i]
def main():
  #initialization
  args=Namespace( perspecies=True, poscar='', write_dos0=False, write_pdos=True, zoomin=None)
  atoms=atom_data()
  dos=dos_data()
  cont=control_knob(args,atoms,dos)
  plt=plot(cont)
  ana=analysis(atoms,dos)

  ##read in overall dos
  io=iodos(cont,plt)
  io.read_poscar()

  x=atoms.positions
  b=atoms.boundary
  b=np.array([b[0,0],b[1,1],b[2,2]])
  bh=b/2.
  #initialize neighborlist
  ngh_list=[]
  for i in range(atoms.natom):
    ngh_list+=[{}]
  #ngh_list=[{}] will lead to a weird problem...

  #build up the neighborlist and second neighbor list
  for i in range(atoms.natom):
    for j in range(i+1,atoms.natom):
      dr=x[j,:]-x[i,:]
      for k in range(3):
        while (abs(dr[k])>bh[k]):
          if (dr[k]>0):
            dr[k]-=b[k]
          else:
            dr[k]+=b[k]
      dist=np.linalg.norm(dr)
      if ( dist < NN):
        ngh_list[i][j]=[dist,dr]
        ngh_list[j][i]=[dist,dr]


  #sort the neighbor and compute Q and G for each atom
  Q=[]
  Q_dagger=[]
  G=np.zeros([atoms.natom,3,3])
  ngh_id=[]
  for i in range(atoms.natom):
    Q+=[[]]
    Q_dagger+=[[]]
    ngh_id+=[[]]
    if ((atoms.constraints[i,0]=='T') or (atoms.constraints[i,0]=='F')):
      if ( atoms.species[i]==1):
        P=[Ceref_n]
      elif ( atoms.species[i]==0):
        P=[O1ref_n,O2ref_n]
      else:
        print "can only handle two types"
        return
      least_missing_id=-1
      least_missing_value=np.max([len(P[ref][:,0]) for ref in range(len(P))])
      angle_store=[]
      Qunsort=[] 

      #first try to find out which pattern matches the best"
      for ref in range(len(P)):
        angle_store+=[{}]
        Qunsort+=[{}] 
        for j in ngh_list[i].keys():
          for k in range(len(P[ref])):
             angle_store[ref][j,k]=np.array(ngh_list[i][j][DX]).dot(P[ref][k,:])/ngh_list[i][j][DIST]
        find_neigh(ngh_list[i],P[ref],angle_store[ref],Qunsort[ref])
        missing_ngh=len(P[ref])-len(Qunsort[ref])
        if (least_missing_value > missing_ngh):
          least_missing_value=missing_ngh
          least_missing_id=ref

      ref=least_missing_id
      #sort the rest,start from the nearest remaining one
      criteria=0.85
      while ((missing_ngh >0) and (criteria>0)):
        find_neigh(ngh_list[i],P[ref],angle_store[ref],Qunsort[ref],criteria,criteria+0.1)
        missing_ngh=len(P[ref])-len(Qunsort[ref])
        criteria-=0.05

      if (missing_ngh > len(ngh_list[i])):
         print "need a larger neighbor list"
         print "not enough neighbor to build Q"
         return
      elif (missing_ngh >0):
         add_nearest_neigh(ngh_list[i],P[ref],angle_store[ref],Qunsort[ref],missing_ngh)

      Q[i]=np.zeros(P[ref].shape)
      ngh_id[i]=[[]]*len(P[ref])
      for k in range(len(P[ref])):
        Q[i][k,:]=Qunsort[ref][k][DX]
        ngh_id[i][k]=Qunsort[ref][k][ID]
      Q_dagger[i],G[i,:,:]=computeG(P[ref],Q[i])
      del angle_store
      del Qunsort
      del P
  del Q

  A=np.zeros([atoms.natom,3,3,3])
  nye_tensor=np.zeros([atoms.natom,3,3])
  #compute the finit different A
  for atomi in range(atoms.natom):
    for i in range(3):
      for m in range(3):
        deltaG_im=np.zeros(len(ngh_id[atomi]))
        for n in range(len(ngh_id[atomi])):
          deltaG_im[n]=G[ngh_id[atomi][n],i,m]
        deltaG_im-=G[atomi,i,m]
        A[atomi,i,m,:]=np.array(Q_dagger[atomi]).dot(deltaG_im)
        for j in range(3):
          for k in range(3):
            nye_tensor[atomi,j,k]-= perm_parity([j,i,m])*A[atomi,i,m,k]
    print atoms.species[atomi],atomi,atoms.positions[atomi,2],' '.join([str(i) for i in nye_tensor[atomi,:].flat])
    
  #  dist=(x[i,0]-17)**2
  #  dist+=(x[i,2]-17)**2
  #  dist=np.sqrt(dist)
  #  if ((ngh_n[i]!=8) and (ngh_n[i]!=4) and (atoms.constraints[i,0]=='T')):
  #     print dist,atoms.species[i],ngh_n[i]
  #io.read_tot_dosfile()
  #ana.bandgap()
  #print ana.peak_weight_center(-32.5,5)

  ##write
  #if cont.write_dos0:
  #    io.write_tot_dosfile()
  #io.delete_tot_dosfile()

  ##read in partial dos
  #if cont.run_pdos:
  #  io.read_pdos()
  #  #plt.pdos()
  #  if cont.write_pdos:
  #      io.write_pdosfile()
  #  io.delete_pdosfile()

def computeG(P,Q):
  Qdagger=np.linalg.pinv(np.array(Q))
  G=np.array(Qdagger.dot(np.array(P)))
  return Qdagger,G

def find_neigh(neighbor,P0,angles,Qunsort,criteria=0.9,criteria_max=1):
  for j in neighbor.keys():
    if ( type(j) is int):
      compare=np.zeros(len(P0))
      for k in range(len(P0)):
        compare[k]=angles[j,k]
      indexk=np.argmax(compare)
      anglemin=compare[indexk]
      if (anglemin>criteria):
         if (indexk not in Qunsort.keys()): #(Q[ref][indexk]==[]):
           Qunsort[indexk]=[neighbor[j][DIST],neighbor[j][DX],anglemin,j]
           del neighbor[j]
         elif ((Qunsort[indexk][DIST] > neighbor[j][DIST]) and (Qunsort[indexk][ANGLE]<criteria_max)):
           neighbor[Qunsort[indexk][ID]] = [Qunsort[indexk][DIST],Qunsort[indexk][DX]]
           Qunsort[indexk]=[neighbor[j][DIST],neighbor[j][DX],anglemin,j]
           del neighbor[j]

def add_nearest_neigh(neighbor,P0,angles,Qunsort,missing_ngh):
  sort_list=sorted(neighbor.iteritems(),key=lambda (k,v):v[0])
  for item in sort_list[:missing_ngh]:
    j=item[0]
    a={}
    for k in range(len(P0)):
       a[j,k]=angles[j,k]
    sort_angle=sorted(a.iteritems(),key=lambda (k,v):-v)
    for angle in sort_angle:
      if ( j in neighbor.keys()):
        indexk=angle[0][1]
        if (indexk not in Qunsort.keys()): #(Q[ref][indexk]==[]):
          Qunsort[indexk]=[neighbor[j][DIST],neighbor[j][DX],angle[1],j]
          del neighbor[j]
  del sort_list

def perm_parity(lst):
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity   

class Namespace:
  def __init__(self, **kwargs):
    self.__dict__.update(kwargs)

if __name__ == '__main__':
  main()

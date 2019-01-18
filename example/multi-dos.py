import os, sys, inspect
cmd_folder = os.path.realpath(os.path.abspath("../"))
if cmd_folder not in sys.path:
  sys.path.insert(0, cmd_folder)
from doslib.control import *
from doslib.iodos import *
from doslib.data import *
from doslib.plot import *
from doslib.analysis import *

import sys
def main():
  #initialization
  for i in range(3):
    argv=[0,'path','example/inputfile/'+str(i+1)+'/',]
    atoms=atom_data()
    dos=dos_data()
    cont=control_knob(argv,atoms,dos)
    plt=plot(cont)
    io=iodos(cont,plt)

    io.read_poscar()
    
    io.read_tot_dosfile()
    ana=analysis(atoms,dos)
    ana.bandgap()
    ana.peak_finder(-32.5,5)

    del atoms,dos,cont,io,plt,ana
    ##write
    #if cont.write_dos0:
    #    io.write_tot_dosfile()
    #io.delete_tot_dosfile()

    ##read in partial dos
    #if cont.run_pdos:
    #  io.read_pdos()
    #  plt.pdos()
    #  if cont.write_pdos:
    #      io.write_pdosfile()
    #  io.delete_pdosfile()

    #debug
    #dumpclean(cont.__dict__)
    #dumpclean(io.__dict__)
    #dumpclean(atoms.__dict__)
    #dumpclean(dos.__dict__)

if __name__ == '__main__':
  main()

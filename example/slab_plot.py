import numpy as np

def plot_atom_start(self):
  rad=12
  scale=1/23.*2.
  atom=self.atom
  zmin=np.amin(atom.positions[:,2])
  zmax=np.amax(atom.positions[:,2])
  lz=zmax-zmin
  zcenter=(zmin+zmax)*0.5
  Oplot=200
  Ceplot=300
  for i in ('zmin', 'zmax', 'lz','zcenter','Oplot','Ceplot','scale'):
     self.plot_atom_dict[i] = locals()[i]

def plot_atom(self,atomi):
  Oplot=self.plot_atom_dict['Oplot']
  Ceplot=self.plot_atom_dict['Ceplot']
  zcenter=self.plot_atom_dict['zcenter']
  zmin=self.plot_atom_dict['zmin']
  scale=self.plot_atom_dict['scale']
  lz=self.plot_atom_dict['lz']

  atom=self.atom
  dos=self.dos

  y=dos.tot[dos.loc_down:dos.loc_up,0]+dos.tot[dos.loc_down:dos.loc_up,1]
  color =self.cm(2.*abs(atom.positions[atomi,2]-zcenter)/lz)  # color will now be an RGBA tuple
  y_off=(atom.positions[atomi,2]-zmin)/lz
  if (atom.species[atomi]==dos.par_element[1]):
    plt.figure(Oplot)
    ymax=np.amax(y)
    pl.plot(dos.Xenergy[dos.loc_down:dos.loc_up],y_off+y/ymax*scale,c=color)
  if (atom.species[atomi]==dos.par_element[3]):
    plt.figure(Ceplot)
    ymax=np.amax(y)
    pl.plot(dos.Xenergy[dos.loc_down:dos.loc_up],y_off+y/ymax*scale,c=color)

def plot_atom_end(self):
  Oplot=self.plot_atom_dict['Oplot']
  Ceplot=self.plot_atom_dict['Ceplot']
  zcenter=self.plot_atom_dict['zcenter']
  zmin=self.plot_atom_dict['zmin']
  scale=self.plot_atom_dict['scale']

  plt.figure(Oplot)
  plt.ylim(0,1+scale)
  plt.show()
  plt.savefig("DOS_O_INSIDE.png")
  plt.close()

  plt.figure(Ceplot)
  plt.ylim(0,1+scale)
  plt.savefig("DOS_Ce_INSIDE.png")
  plt.close()

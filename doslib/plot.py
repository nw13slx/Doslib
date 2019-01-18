import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl

class plot:
  cm = pl.get_cmap('winter')
  plot_atom_dict={}

  plot=plt.plot

  def __init__(self,control):
      self.control=control
      self.atom=control.atom
      self.dos=control.dos
      self.par_symbol=['s','p','d','f']
  def pdos(self):
      if not self.control.run_pdos:
         return 0
      dos=self.dos
      atom=self.atom
      loc_up=dos.loc_up
      loc_down=dos.loc_down
      control=self.control
      simple=self.simple

    #plot each orbital
      if any(x!=-1 for x in dos.par_element):
        par_symbol=self.par_symbol
        y=[]
        line_label=[]
        line_color=['y','y','r','r','b','b','k','k']
        for i in range(4):
          if dos.par_orbital[i] is not None:
            y+=[dos.par_orbital[i][:,0],dos.par_orbital[i][:,1]]
          else:
            y+=[None,None]
          if atom.symbol:
            line_label+=[atom.symbol[dos.par_element[i]]+"_"+par_symbol[i],None]
          else:
            line_label+=["atomtype="+str(dos.par_element[i])+" "+par_symbol[i]+"_orbital",None]
        simple("pdf_orbital",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color

      if (dos.par_element[2] >=0):
        y=[]
        line_label=['d_t2g',None,'d_eg',None]
        line_color=['r','r','b','b']
        y=[dos.d_t2g[:,0],dos.d_t2g[:,1],dos.d_eg[:,0],dos.d_eg[:,1]]
        simple("d_decompose",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color

      if self.control.perspecies:
        y=[]
        line_label=[]
        color_cycle=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
        line_color=[]
        for i in range(atom.ntype):
          if dos.spin:
            y+=[dos.perspecies[:,i,0]]
            line_color+=[color_cycle[i]]
            y+=[dos.perspecies[:,i,1]]
            line_color+=[color_cycle[i]]
            if atom.symbol is not None:
              line_label+=[atom.symbol[i]]
              line_label+=[None]
            else:
              line_label+=["species "+str(i)]
              line_label+=[None]
          else:
            y+=[dos.perspecies[:,i]]
            line_color+=[None]
            if atom.symbol:
              line_label+=[atom.symbol[i]]
            else:
              line_label+=["species "+str(i)]
        simple("perspecies",dos.Xenergy,y,loc_down,loc_up,line_color,line_label)
        del y, line_label, line_color

  def tot_dos(self):
      dos=self.dos
      print dos.spin
      if dos.spin:
        y=[self.dos.dos0[:,0],self.dos.dos0[:,1]]
      else:
        y=[self.dos.dos0]
      line_label=[None,None]
      line_color=['r','b']
      self.simple("DOS-tot",self.dos.Xenergy,y,self.dos.loc_down,self.dos.loc_up,line_color,line_label)
      del y, line_label, line_color

  def simple(self,name,x,y,l_d,l_u,line_color,line_label):
    print "plotting file",name+".png"
    fig=plt.figure(0, figsize=(3.34, 2.5), dpi=300, facecolor='w', edgecolor='k')

    #print "x range: ", x[l_d],x[l_u]
    for i in range(len(y)):
      if y[i] is not None:
        line=plt.plot(x[l_d:l_u],y[i][l_d:l_u],linewidth=1)
        if line_color[i]:
          plt.setp(line,color=line_color[i])
        if line_label[i]:
          plt.setp(line,label=line_label[i])
    if any(x for x in line_label):
      plt.legend(loc=3,fontsize=6)
    if (name=="DOS-tot"):
        self.dos0_extra()
    if (self.control.center_ef==True):
      plt.xlabel("$E-E_\mathrm{VBM}$")
    else:
      plt.xlabel("$E$")
    fig.set_size_inches(3.34,2.5)
    plt.tight_layout()
    fig.savefig(self.control.name+name+".png",dpi=300)
    plt.close()

  def atom_start(self):
      pass
  def atom_during(self,atomi):
      pass
  def atom_end(self):
      pass

  def dos0_extra(self):
    dos=self.dos
    dos0=self.dos.dos0
    ef=self.dos.Xenergy
    fermiN=self.dos.fermiN
    l_d= self.dos.loc_down
    if dos.spin :
      plt.fill_between(ef[l_d:fermiN], 0, dos0[l_d:fermiN,0],facecolor='red')
      plt.fill_between(ef[l_d:fermiN], 0, dos0[l_d:fermiN,1],facecolor='blue')
    else:
      plt.fill_between(ef[l_d:fermiN], 0, dos0[l_d:fermiN],facecolor='red')


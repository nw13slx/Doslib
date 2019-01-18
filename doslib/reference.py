    if (write_PDOS==True):
      write_pdos(nedos,ef,par_orbital,name)
    f.close()
    del par_orbital,d_t2g,d_eg,All,ef,nedos

  #plot the DOS0 file
  plt.figure(0)
  plt.plot(ef,data[:,1],'r-',ef,-data[:,2],'b-')
  if (centerEf==True):
    plt.xlim(-7,7)
    loc_down = np.argmin(abs(ef+7))
    loc_up = np.argmin(abs(ef-7))
    plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
    fermiN = np.argmin(abs(ef))
  else:
    fermiN = np.argmin(abs(ef-efermi))
  if (zoomIn==True):
    plt.xlim(zoomEmin,zoomEmax)
    loc_down = np.argmin(abs(ef-zoomEmin))
    loc_up = np.argmin(abs(ef-zoomEmax))
    plt.ylim(-np.amax(data[loc_down:loc_up,2]), np.amax(data[loc_down:loc_up,1]))
  if (energyshift !=0):
    delta = efermi-energyshift
  else:
    delta = 0
  plt.fill_between((ef+delta)[:fermiN], 0, data[:fermiN,1],facecolor='red')
  plt.fill_between((ef+delta)[:fermiN], 0, -data[:fermiN,2],facecolor='blue')
  pl.show()
  pl.savefig("DOS0.png")
  plt.close()

def write_pdos(nedos,ef,par_orbital,name):
  par_symbol=['s','p','d','f']
  if (symbol!=None):
    for i in range(4):
      par_symbol[i]=symbol[par_element[i]]+"_"+par_symbol[i]
  else:
    for i in range(4):
      par_symbol[i]=str(par_element[i])+"_"+par_symbol[i]

  for i in range(4):
    if (par_element[i]>=0):
      fdos = open(name+par_symbol[i]+"_PDOS.dat", 'w')
      for n in xrange(nedos):
          fdos.write('%10.3e %10.3e %10.3e\n' %(ef[n],par_orbital[i][n,0],par_orbital[i][n,1]))
      fdos.close()

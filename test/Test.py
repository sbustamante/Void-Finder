import numpy as np
import os

catalog = "voids_00"
voids = np.loadtxt( "../%s/void_regions.dat"%(catalog), usecols=(0,1) ).astype(int)
#Lbox
Lbox = 255

#Loading each void
disjoint = []
for i in voids[:,0]:
    if voids[i-1,1] >1:
      indexes = np.loadtxt( "../%s/void_%d.dat"%(catalog,i) )
      #Constructing file for FOF
      os.system("rm fof.tmp")
      fof = open('fof.tmp', 'w')
      #Total number of particles
      fof.write( "%d\n"%(len(indexes)) )
      #Number of DM particles
      fof.write( "%d\n"%(len(indexes)) )
      #Number of Gas particles
      fof.write( "0\n" )
      #Time
      fof.write( "0.1\n" )
      #Nactive
      fof.write( "0\n" )
      #Saving coordinates
      np.savetxt( fof, indexes, fmt="%d\t%d\t%d" )
      fof.close()
      #Running FOF
      os.system("./fofscr/fof -e 2 -m 1 -p %d < ./fof.tmp"%(Lbox) )
      #Detecting disjoint voids
      grps = np.loadtxt( "fof.grp" )
      grps = grps[1:]
      #Storing flag for disjoint voids (1)
      if len(grps)>sum(grps==1):
	  disjoint.append(1)
	  print i, 1
      else:
	  disjoint.append(0)
	  print i, 0
	

os.system("rm fof.tmp")
os.system("rm fof.grp")
disjoint = np.array(disjoint)
#Saving flags of disjoints voids
np.savetxt( "%s.dsj"%(catalog), disjoint, fmt = "%d" )

#List of disjoint voids
voids_id[voids[:,1]]
np.savetxt( "%s.dsj_list"%(catalog),  )
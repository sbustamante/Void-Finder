import numpy as np
from mayavi import mlab as mlab

catalog = "voids_01"
datos = np.loadtxt( "../%s/void_index.dat"%(catalog) )
datos = datos.reshape( (256,256,256) )

def void(i):
    datos2 = np.copy(datos)
    datos2[datos!=i] = 0
    mlab.contour3d(datos2, contours=2)
    #mlab.pipeline.volume(mlab.pipeline.scalar_field(datos), vmin=0)
    mlab.show()
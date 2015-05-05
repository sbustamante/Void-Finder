Void-Finder
===================
**by: Sebastian Bustamante**
**sebastian.bustamante@udea.edu.co**

This code is intended to identidy bulk voids in cosmological simulations. To do so, three different schemes are 
developed and explored here, all of them based upon the watershed transform procedure. The first two schemes are
based on the fractional anisotropy field of two tensorial methods for classifying the cosmic web, the T-web and
the V-web. This leads to the WT-web and the WV-web schemes. The last method is based on classical schemes where 
the density field is taken to be as the underlying tracer field for void identification.

For each scheme is also possible to apply a nth-order median filtering process in order to reduce noise. 
Furthermore, it is possible to activate a boundary removal procedure, where voids with some mean value across 
their common boundaries are merged if such value is above some user-defined threshold (this guarantees a more 
robust hierarchy of voids).

# Elliptic_LCS_2D
### Alireza Hadjighasem (ETH Zurich) 
#### http://www.zfm.ethz.ch/~hadjighasem/index.html

-----------------------------------------------------------------------------
License:

This software is made public for research use only. It may be modified and redistributed under the terms of the GNU General Public License. 

Algorithm:

This code implements theoretical results developed by the Haller Group at ETH Zurich. 
See georgehaller.com or the Wikipedia page on Lagrangian Coherent Structures for more information. 

Citation:

Please cite [1] if you use the code in your own work.

Please also cite the theoretical work underlying my implementation  as follows:

— OrbitDetection.m and eta_tracing.m are implemented based on [2].

— DetectEllipticRegion.m and SettingPoincareSection.m are implemented based on [3]. 

— SingularityDetection.m is implemented based on [5]. 

— cgTensor.m is implemented based on [4].

----------------------------------------------------------------------------- 
References:

[1] A. Hadjighasem, and G. Haller, Geodesic Transport barriers in Jupiter’s Atmosphere: A Video-Based Analysis,  SIAM Review, 58(1): 69-89 (2016). 

[2] G. Haller, FJ. Beron-Vera, Coherent Lagrangian vortices: the black
    holes of turbulence.  J. Fluid Mech. 731 (2013) R4

[3] D. Karrasch, F. Huhn, G. Haller, Automated detection of coherent Lagrangian vortices in two-dimensional unsteady flows. Proc. Royal Society, 471 (2014) 20140639

[4] M. Farazmand, G. Haller, Computing Lagrangian coherent structures from their variational theory. Chaos. 22 (2012) 013128 

[5] M. Farazmand, D. Blazevski, G. Haller, Shearless transport barriers in unsteady two-dimensional flows and maps. Physica D 278–279 (2014) 44–5.

-----------------------------------------------------------------------------

Tested on Matlab R2015b.

Installation notes :

1) After you unzipped the files to mydir, 
   put the Current Directory in Matlab to mydir

2) In the Matlab command prompt,
   type "add_path" to add the necessary folders to the top of the search path

3) You can now try any of the examples

The altimeter products used in this work are produced by SSALTO/DUACS and distributed by AVISO, 
with support from CNES (http://www.aviso.oceanobs.com/duacs). 

Maintained by Alireza Hadjighasem,

alirezah at ethz dot ch

October 30, 2015.

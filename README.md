# injectorCFDModelling
WIP CFD model of injectors using openFOAM 2006  
Currently using interPhaseChangeFoam, from modified cavitating bullet tutorial.  
This uses incompressible fluid model, as seen in CIRA paper.  

Literature:  
https://drive.google.com/drive/u/0/folders/1pLNUjTEDKc12vxf4TFTfJRrxMRk-E3q4  

<h3>notes<h3>
<h4>mesh<h4>
<h5>current mesh (tet4)<h5>  
-refinement near walls, where should be finer towards center as that is where we are concerned with capturing effects.  
-no slip condition should need finer mesh granularity, however i believe nitrous may be near stagnation. Very unsure of this, I'm theorising from UoT's mesh.  
-hexahedral would be preferable. refinement of current tet mesh using snappyHexMesh my produce this result, however this doesnt address the first point.  
    -hexahedral meshes create better results with smaller number of cells I *think*. People who are versed in OF said so.

<h4>solvers<h4>  
<h5>CIRA<h5>  
-interPhaseChangeFoam  
-assumes isothermal & incompressibility  

<h5>UoT<h5>  
-solver not specified, very likely to be reactingMultiPhaseEulerFoam (this is an euler-euler solver with numerical sharpening on the interface)
-they bring to light X. Fan et. al. Model, an analytical model with an empirical correction factor determined from supercrit CO2. very good correlation in non choked flow.  

  
<h5>UoSalento (cavitatingHsolver)<h5>  
-implementation of antoine equations for Psat. This is likely already present in relevant OF solvers.  
-this requires modelling heat transfer, which could increase computational cost drastically.  
-homogenous mixture approach as used also appears to be applicable to injector modelling, as i do not believe there is slip.  
  
  
In my uninformed opinion, I believe the main source of error in CIRA's approach is in their incompressible assumption leading me to believe the eulerian solver approach, even without heat transfer, would have far better accuracy. This is demonstraited in the more accurate results of UoT's paper. 
However, CIRA do appear to capture the choking of the flow. I'm unsure as to whether this is choking occuring from cavitation, or if it's an artificial limit applied. 

# injectorCFDModelling
WIP CFD model of injectors using openFOAM 2006  
Currently using interPhaseChangeFoam, from modified cavitating bullet tutorial.  
This uses incompressible fluid model, as seen in CIRA paper.  

Literature:  
https://drive.google.com/drive/u/0/folders/1pLNUjTEDKc12vxf4TFTfJRrxMRk-E3q4  

<h1>notes</h1>
<h2>mesh</h2>
<h3>current mesh (tet4)</h3>
<p>-refinement near walls, where should be finer towards center as that is where we are concerned with capturing effects.  <br>
-no slip condition should need finer mesh granularity, but seems to be unimportant as has little affect on the validity of the data we're trying to gather. Very unsure of this, I'm theorising from UoT's mesh.  <br>
-hexahedral would be preferable. Refinement of current tet mesh using snappyHexMesh may produce this result, however this doesnt address the first point.  <br>
    -hexahedral meshes create better results with smaller number of cells I <em>think</em>. People who are versed in OF said so.<br></p>

<h3>new mesh (hex/poly)</h3>

<p>-geometry from patches in fusion "ICLR/COTS 10K 2020/Propulsion/cfd" <br>

-following József Nagy's multiphase project, where he creates a mesh from a geometry defined by stl files <br>

-I'm unsure whether it's best to follow this or create the mesh just from a good blockmeshdict, as snappyhexmesh may default to trying to refine around wall boundary conditions<br>

</p>


<h2>solvers</h2>  
<h3>CIRA</h3>  
<p>-interPhaseChangeFoam  <br>
-assumes isothermal & incompressibility  <br></p>

<h3>UoT</h3>  
<p>-solver not specified, very likely to be reactingMultiPhaseEulerFoam (this is an euler-euler solver with numerical sharpening on the interface)<br>
-they bring to light X. Fan et. al. Model, an analytical model with an empirical correction factor determined from supercrit CO2. very good correlation in non choked flow.  <br></p>

<h3>UoSalento (cavitatingHsolver)</h3>  
<p>-implementation of antoine equations for Psat. This is likely already present in relevant OF solvers.  <br>
-this requires modelling heat transfer, which could increase computational cost drastically.  <br>
-homogenous mixture approach as used also appears to be applicable to prediction of mass flow, though tending to underpredict where non homogenous models overpredict.<br>
-inbuilt model considers flow of hydrogen, nitrogen, water. I believe none are particularly applicable just from considering molecular geometry & IMF <br>
-this solver confuses me, as it claims to be isothermal but appears to contains code to recompute the temperature field should the temp fall outside of a predetermined range. <br>

I believe the main source of error in CIRA's approach is in their incompressible assumption leading me to believe the eulerian solver approach, even without heat transfer, would have far better accuracy. This is demonstraited in the more accurate results of UoT's paper. <br>
I am unsure as to why the incompressible assumption appears to lead to underprediction of mass flow. <br>
</p>

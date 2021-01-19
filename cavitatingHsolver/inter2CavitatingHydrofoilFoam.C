/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    interCavitatingHydrofoilFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids with phase-change
    (e.g. cavitation).  Uses a VOF (volume of fluid) phase-fraction based
    interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a
    single momentum equation is solved.

    The set of phase-change models provided are designed to simulate cavitation
    but other mechanisms of phase-change are supported within this solver
    framework.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "CMULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "phaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"
	    //Info<< "prova" << nl << endl;//troubleshooting row
            surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", dimMass/dimTime, 0)
            );

            mixture->correct();

            #include "alphaEqnSubCycle.H"
            interface.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {

	   	/*	**********************************************************************	*/
		/*	In pEqn.H it has been added a switch in order to enable the upper and	*/
		/*	the lower bounds to the p_rgh field. The following expression has been	*/
		/*	adopted: p_rgh = min(max(p_rgh, MinP), MaxP);				*/
		/*	********************************************************************	*/

                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

    	Switch swTeqn_=transportPropertiesDict.lookupOrDefault<Switch>("swTeqn", false);

	if(swTeqn_)
	{
		#include "calcTSatField.H"
		
		const volScalarField& nut_=turbulence->nut();
		turbDT = nut_/turbPr;

		#include "TEqn.H"

		/*	*****************************************************************		*/
		/*	This routine has been added in order to recompute the temperature		*/
		/*	where it goes below a minimum value or above a maximum value.			*/
		/*	The limiting values of T have been tuned by a trial approach.			*/
		/*	*****************************************************************		*/
		T = min(max(T, MinT), MaxT);
		#include "calcPSatField.H"

	    	Switch swSinghal_=transportPropertiesDict.lookupOrDefault<Switch>("swSinghal", false);
		if(swSinghal_)
		{
			const volScalarField& k_=turbulence->k();
			turbP = 0.39 * rho * k_; 
			pvapor = pSat + turbP/2.0;
		} else {
			pvapor = pSat;
		}
		
	} else {
		pvapor = pSat;
	}

        }


	/*	*************************************************************************	*/
	/*	Print on the log file the information about maximum and minimum pressure,	*/
	/*	and the maximum velocity. It is important to monitor the occurrence of   	*/
	/* 	regions characterized by negative pressures that act as triggers for the 	*/
	/*	detachment of the vapour cavity! This produces non-stability of the comp-	*/
	/*	utations!									*/
	/*	*************************************************************************	*/
        Info << "Max pressure: " << max(p).value() << endl; 
        Info << "Min pressure: " << min(p).value() << endl; 
        Info << "Max temperature: " << max(T).value() << endl; 
        Info << "Mix temperature: " << min(T).value() << endl; 
        Info << "Max velocity magnitude: " << max(mag(U)).value() << endl; 
        Info << "Min velocity magnitude: " << min(mag(U)).value() << endl << endl; 

	Info << "Saturation Properties" << endl; 
	Info << "Max TSat: " << max(TSat).value() << endl; 
	Info << "Min TSat: " << min(TSat).value() << endl; 
	Info << "Max pSat: " << max(pSat).value() << endl; 
	Info << "Min pSat: " << min(pSat).value() << endl; 
	Info << "Max pturb from Singhal: " << max(turbP).value() << endl; 
	Info << "Min pturb from Singhal: " << min(turbP).value() << endl; 
	Info << "Max pvapor: " << max(pvapor).value() << endl; 
	Info << "Min pvapor: " << min(pvapor).value() << endl <<endl; 

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "cavitationModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cavitationModel, 0);
    defineRunTimeSelectionTable(cavitationModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitationModel::cavitationModel
(
	const word& type,
    const dictionary& interfaceDict,
    const phaseModel& phase1,
    const phaseModel& phase2
)
:
    interfaceDict_(interfaceDict),
	phase1_(phase1),
	phase2_(phase2),
	pSat_("pSat", dimPressure, interfaceDict_.lookup("pSat")),
	cavitationDict_(interfaceDict_.optionalSubDict(type+"Coeffs"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cavitationModel::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/phase1_.rho() - phase1_*(1.0/phase1_.rho() - 1.0/phase2_.rho()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cavitationModel::vDotP() const
{
    dimensionedScalar pCoeff(1.0/phase1_.rho()  - 1.0/phase2_.rho());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}

bool Foam::cavitationModel::read()
{
	
    cavitationDict_ = interfaceDict_.optionalSubDict(type()+"Coeffs");
    interfaceDict_.lookup("pSat") >> pSat_;

    return true;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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

#include "AntoineConventional.H"
#include "addToRunTimeSelectionTable.H"
#include "function1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(AntoineConventional, 0);
    addToRunTimeSelectionTable(saturationModel, AntoineConventional, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::AntoineConventional::AntoineConventional(const dictionary& dict)
:
    saturationModel(),
    A_("A", dimless, dict),
    B_("B", dimTemperature, dict),
    C_("C", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::saturationModels::AntoineConventional::~AntoineConventional()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineConventional::pSat
(
    const volScalarField& T
) const
{
    return
        dimensionedScalar("one", dimPressure, 1)
       *100000*pow(10,(A_ - B_/(C_ + T)));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineConventional::pSatPrime
(
    const volScalarField& T
) const
{
    return  pSat(T)*B_*log(dimensionedScalar("one", dimless, 1)*10)
	   /pow((C_ + T),2);
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineConventional::lnPSat
(
    const volScalarField& T
) const
{
    return log(pSat(T));
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::AntoineConventional::Tsat
(
    const volScalarField& p
) const
{
    return
        B_/(A_ - log10(p*dimensionedScalar("one", dimless/dimPressure, 1)))
      - C_;
}


// ************************************************************************* //

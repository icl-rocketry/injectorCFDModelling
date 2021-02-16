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

#include "SchnerrSauer.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace cavitationModels
{
    defineTypeNameAndDebug(SchnerrSauer, 0);
    addToRunTimeSelectionTable
    (
        cavitationModel,
        SchnerrSauer,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cavitationModels::SchnerrSauer::SchnerrSauer
(
    const dictionary& interfaceDict,
    const phaseModel& phase1,
	const phaseModel& phase2
)
:
    cavitationModel(typeName, interfaceDict, phase1, phase2),

    n_("n", dimless/dimVolume, cavitationDict_),
    dNuc_("dNuc", dimLength, cavitationDict_),
    Cc_("Cc", dimless, cavitationDict_),
    Cv_("Cv", dimless, cavitationDict_),

    p0_("0", pSat().dimensions(), 0.0)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::cavitationModels::SchnerrSauer::rRb
(
    const volScalarField& limitedAlpha1
) const
{
    return pow
    (
        ((4*constant::mathematical::pi*n_)/3)
       *limitedAlpha1/(1.0 + alphaNuc() - limitedAlpha1),
        1.0/3.0
    );
}


Foam::dimensionedScalar
Foam::cavitationModels::SchnerrSauer::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}


Foam::tmp<Foam::volScalarField>
Foam::cavitationModels::SchnerrSauer::pCoeff
(
    const volScalarField& p
) const
{
    volScalarField limitedAlpha1(min(max(phase1_, scalar(0)), scalar(1)));
    volScalarField rho
    (
        limitedAlpha1*phase1_.rho() + (scalar(1) - limitedAlpha1)*phase2_.rho()
    );

    return
        (3*phase1_.rho()*phase2_.rho())*sqrt(2/(3*phase1_.rho()))
       *rRb(limitedAlpha1)/(rho*sqrt(mag(p - pSat()) + 0.01*pSat()));
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cavitationModels::SchnerrSauer::mDotAlphal() const
{
    const volScalarField& p = phase1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));
    volScalarField limitedAlpha1(min(max(phase1_, scalar(0)), scalar(1)));

    return Pair<tmp<volScalarField>>
    (
        Cc_*limitedAlpha1*pCoeff*max(p - pSat(), p0_),

        Cv_*(1.0 + alphaNuc() - limitedAlpha1)*pCoeff*min(p - pSat(), p0_)
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::cavitationModels::SchnerrSauer::mDotP() const
{
    const volScalarField& p = phase1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(phase1_, scalar(0)), scalar(1)));
    volScalarField apCoeff(limitedAlpha1*pCoeff);

    return Pair<tmp<volScalarField>>
    (
        Cc_*(1.0 - limitedAlpha1)*pos0(p - pSat())*apCoeff,

        (-Cv_)*(1.0 + alphaNuc() - limitedAlpha1)*neg(p - pSat())*apCoeff
    );
}

bool Foam::cavitationModels::SchnerrSauer::read()
{
    if (cavitationModel::read())
    {
        cavitationDict_ = interfaceDict_.optionalSubDict(type()+"Coeffs");

        cavitationDict_.lookup("n") >> n_;
        cavitationDict_.lookup("dNuc") >> dNuc_;
        cavitationDict_.lookup("Cc") >> Cc_;
        cavitationDict_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //

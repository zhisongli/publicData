/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "kOmegaSSTNew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// -------------------------------------------------------------------------
// correctNut — 3-argument stabilized version (NEW)
//
// This overload is called by kOmegaSSTNewBase::correct() when the
// Larsen & Fuhrman (2018) stabilization is active.  The wrapper
// delegates to the base class to set nut_, then calls
// BasicTurbulenceModel::correctNut() to update the turbulence
// thermal diffusivity (alphat) — exactly the same pattern as the
// original 1-arg overload below.
// -------------------------------------------------------------------------
template<class BasicTurbulenceModel>
void kOmegaSSTNew<BasicTurbulenceModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& pOmega,
    const volScalarField& F23
)
{
    // Correct the turbulence viscosity (stabilized limiter in base class)
    kOmegaSSTNewBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2,
        pOmega,
        F23
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


// -------------------------------------------------------------------------
// correctNut — 1-argument standard version (UNCHANGED from v2512)
//
// Called by kOmegaSSTNewBase::correct() when stabilization is OFF.
// -------------------------------------------------------------------------
template<class BasicTurbulenceModel>
void kOmegaSSTNew<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    // Correct the turbulence viscosity
    kOmegaSSTNewBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>::correctNut
    (
        S2
    );

    // Correct the turbulence thermal diffusivity
    BasicTurbulenceModel::correctNut();
}


// -------------------------------------------------------------------------
// correctNut — zero-argument fallback (UNCHANGED from v2512)
//
// Computes S2 from gradU internally and delegates to the 1-arg version.
// Note: this does NOT need to handle the stabilization routing itself —
// the base-class zero-arg correctNut() already dispatches to the right
// overload.  But we keep this override so that the thermal diffusivity
// correction is guaranteed to happen regardless of the call path.
// -------------------------------------------------------------------------
template<class BasicTurbulenceModel>
void kOmegaSSTNew<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTNew<BasicTurbulenceModel>::kOmegaSSTNew
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSSTNewBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //

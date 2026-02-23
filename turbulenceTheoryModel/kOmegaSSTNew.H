/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::RASModels::kOmegaSSTNew

Group
    grpRASTurbulence

Description
    Implementation of the k-omega-SST turbulence model for
    incompressible and compressible flows.

    Light wrapper around base class.

    Turbulence model described in:
    \verbatim
        Menter, F. R. & Esch, T. (2001).
        Elements of Industrial Heat Transfer Prediction.
        16th Brazilian Congress of Mechanical Engineering (COBEM).
    \endverbatim

    with updated coefficients from
    \verbatim
        Menter, F. R., Kuntz, M., and Langtry, R. (2003).
        Ten Years of Industrial Experience with the SST Turbulence Model.
        Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano,
        & M. Tummers, Begell House, Inc., 625 - 632.
    \endverbatim

    but with the consistent production terms from the 2001 paper as form in the
    2003 paper is a typo, see
    \verbatim
        http://turbmodels.larc.nasa.gov/sst.html
    \endverbatim

    and the addition of the optional F3 term for rough walls from
    \verbatim
        Hellsten, A. (1998).
        "Some Improvements in Menter's k-omega-SST turbulence model"
        29th AIAA Fluid Dynamics Conference, AIAA-98-2554.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficuients in a consistent
    manner.  The paper suggests that sigma is blended but this would not be
    consistent with the blending of the k-epsilon and k-omega models.

    Also note that the error in the last term of equation (2) relating to
    sigma has been corrected.

    Wall-functions are applied in this implementation by using equations (14)
    to specify the near-wall omega as appropriate.

    The blending functions (15) and (16) are not currently used because of the
    uncertainty in their origin, range of applicability and that if y+ becomes
    sufficiently small blending u_tau in this manner clearly becomes nonsense.

    The default model coefficients are
    \verbatim
        kOmegaSSTNewBaseCoeffs
        {
            alphaK1         0.85;
            alphaK2         1.0;
            alphaOmega1     0.5;
            alphaOmega2     0.856;
            beta1           0.075;
            beta2           0.0828;
            betaStar        0.09;
            gamma1          5/9;
            gamma2          0.44;
            a1              0.31;
            b1              1.0;
            c1              10.0;
            F3              no;

            // Optional decay control (Spalart & Rumsey 2007)
            decayControl    yes;
            kInf            \<far-field k value\>;
            omegaInf        \<far-field omega value\>;

            // Optional Larsen & Fuhrman (2018) stabilization
            stabilization   yes;
            lambda2         0.05;   // Stabilization blending coefficient
            alphaBS         1.36;   // Buoyancy-production multiplier

            // Optional Hellsten (1998) Rotation/Curvature Correction
            controlRCHellsten   no;
            cRC             1.4;

            // Optional Spalart-Shur / Menter-Smirnov (2009) Curvature Correction
            // (mutually exclusive with controlRCHellsten)
            controlRCMentor     no;
            cr1             1.0;
            cr2             2.0;
            cr3             1.0;
        }
    \endverbatim

SourceFiles
    kOmegaSSTNew.C

SeeAlso
    kOmegaSSTNewBase.H

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTNew_H
#define kOmegaSSTNew_H

#include "kOmegaSSTNewBase.H"
#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                          Class kOmegaSSTNew Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTNew
:
    public kOmegaSSTNewBase<eddyViscosity<RASModel<BasicTurbulenceModel>>>
{
    // Private Member Functions

        //- No copy construct
        kOmegaSSTNew(const kOmegaSSTNew&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTNew&) = delete;


protected:

    // Protected Member Functions

        // -----------------------------------------------------------------
        // Three correctNut overloads mirror the base class.
        //
        // The wrapper's ONLY job is to additionally call
        // BasicTurbulenceModel::correctNut() after the base-class
        // implementation has set nut_, so that the turbulence thermal
        // diffusivity (alphat) is also updated.
        //
        // (1) correctNut(S2, pOmega, F23) — stabilized path (new)
        // (2) correctNut(S2)              — standard path   (original)
        // (3) correctNut()                — zero-arg fallback
        // -----------------------------------------------------------------

        //- Stabilized nut correction (Larsen & Fuhrman 2018)
        virtual void correctNut
        (
            const volScalarField& S2,
            const volScalarField& pOmega,
            const volScalarField& F23
        );

        //- Standard nut correction from S2
        virtual void correctNut(const volScalarField& S2);

        //- Zero-argument fallback
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTNew");


    // Constructors

        //- Construct from components
        kOmegaSSTNew
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTNew() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
    #include "kOmegaSSTNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //

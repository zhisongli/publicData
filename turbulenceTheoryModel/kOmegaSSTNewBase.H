/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::kOmegaSSTNewBase

Description
    Base class implementation of the k-omega-SST turbulence model for
    incompressible and compressible flows, with four optional extensions.

    This model is directly implemented from the OpenFOAM-v2512. The viscosity
    limitter, the buoyant effect and the rotation correction are added.

    -----------------------------------------------------------------------
    EXTENSION 1 — Turbulence viscosity limitter
    -----------------------------------------------------------------------
    Adds a third denominator candidate in the eddy-viscosity expression
    that scales with the ratio S2/pOmega (strain-rate squared / vorticity
    squared).  In wave-orbital-motion regions where S ≈ Ω the limiter
    activates and suppresses nut, preventing the well-known spurious
    turbulence growth beneath surface waves.  In shear-dominated boundary
    layers S2/pOmega ≈ 1 and the limiter is dormant.

    Activated by setting "stabilization yes"; in the turbulence dictionary.

    \verbatim
        Larsen, B. E. and Fuhrman, D. R. (2018).
        On the over-production of turbulence beneath surface waves
        in Reynolds-averaged Navier-Stokes models.
        J. Fluid Mech., 853, 419-460.
    \endverbatim

    -----------------------------------------------------------------------
    EXTENSION 2 — Buoyancy production term
    -----------------------------------------------------------------------
    Automatically detected and activated at run-time whenever both the
    mixture density field "rho" and the gravity vector "g" exist in the
    objectRegistry (i.e. in multiphase VOF solvers).  Zero overhead in
    single-phase runs.

    Under stable stratification (water below air):
      g · ∇ρ > 0  →  implicit sink on k  →  suppresses turbulence  ✓
    Under unstable stratification:
      g · ∇ρ < 0  →  explicit source on k  →  enhances mixing  ✓

    \verbatim
        Devolder, B., Rauwoens, P., & Troch, P. (2017). 
        Application of a buoyancy-modified k-ω SST turbulence model 
        to simulate wave run-up around a monopile subjected to 
        regular waves using OpenFOAM®. 
        Coastal Engineering, 125, 81-94.

        Devolder, B., Troch, P., & Rauwoens, P. (2018). 
        Performance of a buoyancy-modified k-ω and k-ω SST turbulence 
        model for simulating wave breaking under regular waves using 
        OpenFOAM®. 
        Coastal Engineering, 138, 49-65.

        Larsen, B. E. and Fuhrman, D. R. (2018).
        On the over-production of turbulence beneath surface waves
        in Reynolds-averaged Navier-Stokes models.
        J. Fluid Mech., 853, 419-460.        
    \endverbatim

    -----------------------------------------------------------------------
    EXTENSION 3 — Rotation/Curvature Correction (RCH) to the omega equation
    -----------------------------------------------------------------------
    Modifies the destruction term in the omega equation by multiplying
    it with a correction factor F4 = 1/(1 + cRC * Ri), where Ri is
    the Richardson number derived from the ratio of rotation rate to
    strain rate:  Ri = (W/S) * (W/S - 1).

    When flow is curved or rotating, W ≠ S and the Richardson number
    becomes nonzero:
      Ri > 0 (stabilising rotation, W > S)  →  F4 < 1  →  less ω
        destruction  →  more ω  →  less k  →  suppresses turbulence  ✓
      Ri < 0 (destabilising, W < S)  →  F4 > 1  →  more ω destruction
        →  less ω  →  more k  →  enhances turbulence  ✓

    Activated by setting "controlRCHellsten yes;".

    \verbatim
        Hellsten, A. (1998).
        Some Improvements in Menter''s k-omega SST Turbulence Model.
        AIAA Paper 98-2554.

        Mani, M., Ladd, J. A. and Bower, W. W. (2004).
        Rotation and Curvature Correction Assessment for One- and
        Two-Equation Turbulence Models.
        Journal of Aircraft, 41(2), 268-273.
    \endverbatim

    -----------------------------------------------------------------------
    EXTENSION 4 — Spalart-Shur / Menter-Smirnov (2009) Curvature
                   Correction (CC)
    -----------------------------------------------------------------------
    Modifies the production term in both k and omega equations via a
    multiplier Fr1 that depends on the Lagrangian (material) derivative
    of the strain-rate tensor.  This captures streamline curvature and
    system rotation effects more generally than the Hellsten approach.

    Fr1 = max( min( f_rotation, 1.25 ), 0 )
    where f_rotation = (1+cr1) * 2r*\/(1+r*) * (1 - cr3*atan(cr2*r̂)) - cr1
      r*    = |S| / |W|          (strain-to-rotation ratio)
      r̂     = 2*W_ik*S_jk*(DS_ij/Dt) / (|W| * D³)
      D     = sqrt(max(S², 0.09ω²))

    Activated by setting "controlRCMentor yes;".

    NOTE: controlRCHellsten and controlRCMentor are MUTUALLY EXCLUSIVE.
    Enabling both will cause a FatalError at construction time.

    \verbatim
        Smirnov, P. E. and Menter, F. R. (2009).
        Sensitization of the SST Turbulence Model to Rotation and
        Curvature by Applying the Spalart-Shur Correction Term.
        ASME J. Turbomachinery, 131, 041010.
    \endverbatim

    -----------------------------------------------------------------------

    Base SST model described in:
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

    and the optional decay control from:
    \verbatim
        Spalart, P. R. and Rumsey, C. L. (2007).
        Effective Inflow Conditions for Turbulence Models in Aerodynamic
        Calculations.
        AIAA Journal, 45(10), 2544 - 2553.
    \endverbatim

    Note that this implementation is written in terms of alpha diffusion
    coefficients rather than the more traditional sigma (alpha = 1/sigma) so
    that the blending can be applied to all coefficients in a consistent
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

Author
    Larsen and Fuhrman, Technical University of Denmark (DTU)
    Zhisong Li, Shanghai Jiao Tong University (SJTU) and
    Shanghai Maritime University (SMU)

SourceFiles
    kOmegaSSTNewBase.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTNewBase_H
#define kOmegaSSTNewBase_H

// Required for dynamic look-up of the gravity vector "g" from the registry
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTNewBase Declaration
\*---------------------------------------------------------------------------*/

template<class BasicEddyViscosityModel>
class kOmegaSSTNewBase
:
    public BasicEddyViscosityModel
{
protected:

    // Protected Data

        // Standard SST model coefficients

            dimensionedScalar alphaK1_;
            dimensionedScalar alphaK2_;

            dimensionedScalar alphaOmega1_;
            dimensionedScalar alphaOmega2_;

            dimensionedScalar gamma1_;
            dimensionedScalar gamma2_;

            dimensionedScalar beta1_;
            dimensionedScalar beta2_;

            dimensionedScalar betaStar_;

            dimensionedScalar a1_;
            dimensionedScalar b1_;
            dimensionedScalar c1_;

            //- Flag to include the F3 rough-wall correction (Hellsten 1998)
            Switch F3_;


        // Fields

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;

            //- Turbulent kinetic energy field [m^2/s^2]
            volScalarField k_;

            //- Specific dissipation rate field [1/s]
            volScalarField omega_;


        // Decay control (Spalart & Rumsey 2007)

            //- Flag to include the decay control
            Switch decayControl_;
            dimensionedScalar kInf_;
            dimensionedScalar omegaInf_;


        // ------------------------------------------------------------------
        // Larsen & Fuhrman (2018) stabilization coefficients
        //
        //   stabilization_  : master on/off switch (default: off, so the
        //                     standard SST model is recovered by default)
        //   lambda2_        : weighting factor (eq. 21 in L&F 2018).
        //                     Controls how aggressively the strain-rate
        //                     limiter acts.  Default 0.05.
        //   alphaBS_        : multiplier for the buoyancy production term
        //                     Gb = -alphaBS * nut/sigmaT * (g . grad(rho)).
        //                     Default 1.36 (Devolder et al. 2017, 2018).
        //   pOmegaSmall_    : small number [1/s^2] added to the denominator
        //                     of the stabilization limiter to avoid division
        //                     by zero in irrotational regions.
        // ------------------------------------------------------------------

            Switch stabilization_;
            dimensionedScalar lambda2_;
            dimensionedScalar alphaBS_;
            dimensionedScalar pOmegaSmall_;


        // ------------------------------------------------------------------
        // Hellsten (1998) Rotation/Curvature Correction (RCH)
        //
        //   controlRCHellsten_ : on/off switch (default: off)
        //   cRC_               : Richardson-number sensitivity coefficient
        //                        (default 1.4, Mani et al. 2004).
        //
        //   When active, a correction factor F4 = 1/(1 + cRC * Ri) is
        //   applied to the omega destruction term, where the Richardson
        //   number is Ri = (W/S)*(W/S - 1).
        //
        //   - Stabilising rotation (W > S): Ri > 0 → F4 < 1
        //     → ω destruction reduced → ω increases → k suppressed  ✓
        //   - Destabilising curvature (W < S): Ri < 0 → F4 > 1
        //     → ω destruction enhanced → ω decreases → k boosted  ✓
        // ------------------------------------------------------------------

            Switch controlRCHellsten_;
            dimensionedScalar cRC_;


        // ------------------------------------------------------------------
        // Spalart-Shur / Menter-Smirnov (2009) Curvature Correction (CC)
        //
        //   controlRCMentor_ : on/off switch (default: off)
        //   cr1_, cr2_, cr3_ : model constants (defaults 1.0, 2.0, 1.0)
        //
        //   When active, a multiplier Fr1 is applied to the production
        //   terms in BOTH the k and omega equations:
        //     Pk  → Pk  * Fr1
        //     Pω  → Pω  * Fr1
        //
        //   Fr1 is computed from the Lagrangian derivative of the
        //   strain-rate tensor, capturing streamline curvature and
        //   system rotation effects.
        //
        //   NOTE: Mutually exclusive with controlRCHellsten.
        // ------------------------------------------------------------------

            Switch controlRCMentor_;
            dimensionedScalar cr1_;
            dimensionedScalar cr2_;
            dimensionedScalar cr3_;


    // Protected Member Functions

        //- Set decay control with kInf and omegaInf
        void setDecayControl(const dictionary& dict);

        virtual tmp<volScalarField> F1
        (
            const volScalarField& CDkOmega
        ) const;
        virtual tmp<volScalarField> F2() const;
        virtual tmp<volScalarField> F3() const;
        virtual tmp<volScalarField> F23() const;

        //- Return the blended field
        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        //- Return the internal blended field
        tmp<volScalarField::Internal> blend
        (
            const volScalarField::Internal& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
            return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaK
        (
            const volScalarField& F1
        ) const
        {
            return blend
            (
                F1, 
                alphaK1_, 
                alphaK2_
            );
        }

        tmp<volScalarField> alphaOmega
        (
            const volScalarField& F1
        ) const
        {
            return blend
            (
                F1, 
                alphaOmega1_, 
                alphaOmega2_
            );
        }

        tmp<volScalarField::Internal> beta
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                IOobject::scopedName
                (
                    this->type(), 
                    "beta"
                ),
                blend(F1, beta1_, beta2_)
            );
        }

        tmp<volScalarField::Internal> gamma
        (
            const volScalarField::Internal& F1
        ) const
        {
            return tmp<volScalarField::Internal>::New
            (
                IOobject::scopedName
                (
                    this->type(), 
                    "gamma"
                ),
                blend(F1, gamma1_, gamma2_)
            );
        }

        // -----------------------------------------------------------------
        // correctNut overloads
        //
        // (1) Three-argument: stabilized path (S2, pOmega, F23).
        //     Primary workhorse when stabilization is ON.
        //
        // (2) One-argument: standard path (S2 only).
        //     Used when stabilization is OFF.
        //
        // (3) Zero-argument: fallback that computes gradU internally
        //     and routes to the appropriate overload.
        // -----------------------------------------------------------------

        virtual void correctNut // new added by ZhisongLi
        (
            const volScalarField& S2,
            const volScalarField& pOmega,
            const volScalarField& F23
        );

        virtual void correctNut
        (
            const volScalarField& S2
        );

        virtual void correctNut();

        //- Return square of strain rate
        virtual tmp<volScalarField> S2
        (
            const volTensorField& gradU
        ) const;

        //- Return k production rate
        virtual tmp<volScalarField::Internal> Pk
        (
            const volScalarField::Internal& G
        ) const;

        //- Return epsilon/k which for standard RAS is betaStar*omega
        virtual tmp<volScalarField::Internal> epsilonByk
        (
            const volScalarField& F1,
            const volTensorField& gradU
        ) const;

        //- Return (G/nu)_0
        virtual tmp<volScalarField::Internal> GbyNu0
        (
            const volTensorField& gradU,
            const volScalarField& S2
        ) const;

        //- Return G/nu
        virtual tmp<volScalarField::Internal> GbyNu
        (
            const volScalarField::Internal& GbyNu0,
            const volScalarField::Internal& F2,
            const volScalarField::Internal& S2
        ) const;

        virtual tmp<fvScalarMatrix> kSource() const;

        virtual tmp<fvScalarMatrix> omegaSource() const;

        virtual tmp<fvScalarMatrix> Qsas
        (
            const volScalarField::Internal& S2,
            const volScalarField::Internal& gamma,
            const volScalarField::Internal& beta
        ) const;


public:

    typedef typename BasicEddyViscosityModel::alphaField alphaField;
    typedef typename BasicEddyViscosityModel::rhoField rhoField;
    typedef typename BasicEddyViscosityModel::transportModel transportModel;


    // Generated Methods

        //- No copy construct
        kOmegaSSTNewBase(const kOmegaSSTNewBase&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTNewBase&) = delete;


    // Constructors

        //- Construct from components
        kOmegaSSTNewBase
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~kOmegaSSTNewBase() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DkEff",
                alphaK(F1)*this->nut_ + this->nu()
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const
        {
            return tmp<volScalarField>::New
            (
                "DomegaEff",
                alphaOmega(F1)*this->nut_ + this->nu()
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTNewBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

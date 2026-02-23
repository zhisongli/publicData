/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "kOmegaSSTNewBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// =========================================================================
// Blending functions F1, F2, F3, F23
// (Unchanged from standard SST — Menter 2001/2003)
// =========================================================================

template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTNewBase<BasicEddyViscosityModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    // CDkOmegaPlus: clamp cross-diffusion from below to prevent tanh(0/0).
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar(dimless/sqr(dimTime), 1.0e-10)
    );

    // arg1 = min(min(max(√k/(β*ω*y), 500ν/(ω*y²)), 4αω2*k/(CD*y²)), 10)
    // F1 = tanh(arg1⁴): equals 1 near walls (k-ω zone), 0 in free stream
    // (k-ε zone).
    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTNewBase<BasicEddyViscosityModel>::F2() const
{
    // F2 = tanh(arg2²): equals 1 in boundary layer, 0 in free stream.
    // Used in the Bradshaw limiter for nut.
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*(this->mu()/this->rho_)/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTNewBase<BasicEddyViscosityModel>::F3() const
{
    // Optional rough-wall correction (Hellsten 1998).
    // Only activated when the "F3" switch is true.
    tmp<volScalarField> arg3 = min
    (
        150*(this->mu()/this->rho_)/(omega_*sqr(y_)),
        scalar(10)
    );

    return 1 - tanh(pow4(arg3));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField> kOmegaSSTNewBase<BasicEddyViscosityModel>::F23() const
{
    tmp<volScalarField> f23(F2());

    if (F3_)
    {
        f23.ref() *= F3();
    }

    return f23;
}


// =========================================================================
// correctNut — 3-argument stabilized version (Larsen & Fuhrman 2018)
//
//   Standard SST:
//     nut = a1 * k / max(a1*omega, b1*F23*sqrt(S2))
//
//   Stabilized SST (eq. 21 in L&F 2018):
//     nut = a1 * k / max( a1 * lambda2*(beta1/(betaStar*gamma1))
//                              * S2/(pOmega+small) * omega,
//                         max(a1*omega, b1*F23*sqrt(S2)) )
//
//   The L&F limiter term activates only when S2/pOmega >> 1 (irrotational
//   wave regions).  In boundary layers S2/pOmega ≈ 1 so the Bradshaw
//   limiter remains dominant and wall-bounded flow is unaffected.
//
//   pOmega = 2*|W|² (W = skew(gradU), the rotation-rate tensor).
// =========================================================================
template<class BasicEddyViscosityModel>
void kOmegaSSTNewBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2,
    const volScalarField& pOmega,
    const volScalarField& F23
)
{
    this->nut_ = a1_*k_/max
    (
        // --- Larsen & Fuhrman rotation-rate limiter ---
        a1_*lambda2_*beta1_/(betaStar_*gamma1_)
            * S2 / max(pOmega, pOmegaSmall_) * omega_, // more reasonable
            // * S2 / (pOmega + pOmegaSmall_) * omega_,

        // --- Standard Bradshaw limiter ---
        max
        (
            a1_*omega_, 
            b1_*F23*sqrt(S2)
        )
    );

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


// =========================================================================
// correctNut — 1-argument version (standard path when stabilization OFF)
//
//   Identical to the original v2512 SST formula.  Retained for backward
//   compatibility and used as the default path.
// =========================================================================
template<class BasicEddyViscosityModel>
void kOmegaSSTNewBase<BasicEddyViscosityModel>::correctNut
(
    const volScalarField& S2
)
{
    // Correct the turbulence viscosity
    this->nut_ = a1_*k_/max(a1_*omega_, b1_*F23()*sqrt(S2));
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
}


// =========================================================================
// correctNut — zero-argument fallback
//
//   Called when no gradU information is available (e.g. base-class init).
//   Routes to the appropriate overload.
// =========================================================================
template<class BasicEddyViscosityModel>
void kOmegaSSTNewBase<BasicEddyViscosityModel>::correctNut()
{
    if (stabilization_)
    {
        tmp<volTensorField> tgradU = fvc::grad(this->U_);
        tmp<volScalarField> tS2     = 2.0*magSqr(symm(tgradU()));
        tmp<volScalarField> tpOmega = 2.0*magSqr(skew(tgradU()));
        correctNut(tS2(), tpOmega(), F23());
    }
    else
    {
        correctNut(2*magSqr(symm(fvc::grad(this->U_))));
    }
}


// =========================================================================
// Helper functions: S2, GbyNu0, GbyNu, Pk, epsilonByk
// (Unchanged from standard SST)
// =========================================================================

template<class BasicEddyViscosityModel>
Foam::tmp<Foam::volScalarField> 
kOmegaSSTNewBase<BasicEddyViscosityModel>::S2
(
    const volTensorField& gradU
) const
{
    // S2 = 2*S_ij*S_ij, where S = symm(gradU)
    return 2*magSqr(symm(gradU));
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal> 
kOmegaSSTNewBase<BasicEddyViscosityModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    // Production limiter: min(G, c1*betaStar*k*omega)
    // Prevents runaway production in high-shear / stagnation zones.
    return min(G, (c1_*betaStar_)*this->k_()*this->omega_());
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTNewBase<BasicEddyViscosityModel>::epsilonByk
(
    const volScalarField& /* F1 not used */,
    const volTensorField& /* gradU not used */
) const
{
    // Standard SST specific dissipation: epsilon/k = betaStar * omega
    return betaStar_*omega_();
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTNewBase<BasicEddyViscosityModel>::GbyNu0
(
    const volTensorField& gradU,
    const volScalarField& /* S2 not used */
) const
{
    // Un-limited specific production: G/nu = gradU:(dev(2*symm(gradU)))
    return tmp<volScalarField::Internal>::New
    (
        IOobject::scopedName
        (
            this->type(), 
            "GbyNu"
        ),
        gradU() && devTwoSymm(gradU())
    );
}


template<class BasicEddyViscosityModel>
tmp<volScalarField::Internal>
kOmegaSSTNewBase<BasicEddyViscosityModel>::GbyNu
(
    const volScalarField::Internal& GbyNu0,
    const volScalarField::Internal& F2,
    const volScalarField::Internal& S2
) const
{
    // Bradshaw / stagnation-point limiter.
    // Clips G/nu to c1/a1 * betaStar * omega * max(a1*omega, b1*F2*sqrt(S2)).
    // Prevents over-production in impingement and stagnation regions.
    return min
    (
        GbyNu0,
        (c1_/a1_)*betaStar_*omega_()*max(a1_*omega_(), b1_*F2*sqrt(S2))
    );
}


// =========================================================================
// Source term stubs (zero by default; override in derived models e.g. SAS)
// =========================================================================

template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix>
kOmegaSSTNewBase<BasicEddyViscosityModel>::kSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        k_,
        dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix>
kOmegaSSTNewBase<BasicEddyViscosityModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}


template<class BasicEddyViscosityModel>
tmp<fvScalarMatrix> kOmegaSSTNewBase<BasicEddyViscosityModel>::Qsas
(
    const volScalarField::Internal& S2,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    return tmp<fvScalarMatrix>::New
    (
        omega_,
        dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
kOmegaSSTNewBase<BasicEddyViscosityModel>::kOmegaSSTNewBase
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicEddyViscosityModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    // --- Standard SST coefficients ---
    alphaK1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK1",
            this->coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaK2",
            this->coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega1",
            this->coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaOmega2",
            this->coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma1",
            this->coeffDict_,
            5.0/9.0
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "gamma2",
            this->coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta1",
            this->coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "beta2",
            this->coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "a1",
            this->coeffDict_,
            0.31
        )
    ),
    b1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "b1",
            this->coeffDict_,
            1.0
        )
    ),
    c1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "c1",
            this->coeffDict_,
            10.0 // in NASA, c1 = 20.
        )
    ),
    F3_
    (
        Switch::getOrAddToDict
        (
            "F3",
            this->coeffDict_,
            false
        )
    ),

    y_(wallDist::New(this->mesh_).y()),

    k_
    (
        IOobject
        (
            IOobject::groupName
            (
                "k", 
                alphaRhoPhi.group()
            ),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName
            (
                "omega", 
                alphaRhoPhi.group()
            ),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    // --- Decay control (Spalart & Rumsey 2007) ---
    decayControl_
    (
        Switch::getOrAddToDict
        (
            "decayControl",
            this->coeffDict_,
            false
        )
    ),
    kInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kInf",
            this->coeffDict_,
            k_.dimensions(),
            0
        )
    ),
    omegaInf_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaInf",
            this->coeffDict_,
            omega_.dimensions(),
            0
        )
    ),

    // --- Larsen & Fuhrman (2018) stabilization ---
    stabilization_
    (
        Switch::getOrAddToDict
        (
            "stabilization",
            this->coeffDict_,
            false           // OFF by default — standard SST behaviour
        )
    ),
    lambda2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "lambda2",
            this->coeffDict_,
            0.05            // Default from Larsen & Fuhrman (2018)
        )
    ),
    alphaBS_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaBS",
            this->coeffDict_,
            1.36            // Default from Larsen & Fuhrman (2018)
        )
    ),
    // pOmegaSmall_("pOmegaSmall", dimless/sqr(dimTime), SMALL),
    pOmegaSmall_("pOmegaSmall", dimless/sqr(dimTime), 1e-6), // less like divergence

    // --- Hellsten (1998) Rotation/Curvature Correction ---
    controlRCHellsten_
    (
        Switch::getOrAddToDict
        (
            "controlRCHellsten",
            this->coeffDict_,
            false
        )
    ),
    cRC_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cRC",
            this->coeffDict_,
            1.4             // Default from Mani et al. (2004)
        )
    ),

    // --- Spalart-Shur / Menter-Smirnov (2009) Curvature Correction ---
    controlRCMentor_
    (
        Switch::getOrAddToDict
        (
            "controlRCMentor",
            this->coeffDict_,
            false
        )
    ),
    cr1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr1",
            this->coeffDict_,
            1.0
        )
    ),
    cr2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr2",
            this->coeffDict_,
            2.0
        )
    ),
    cr3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "cr3",
            this->coeffDict_,
            1.0
        )
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    setDecayControl(this->coeffDict_);

    // -----------------------------------------------------------------
    // Mutual exclusivity check: Hellsten RCH and Menter-Smirnov CC
    // cannot both be active because they modify different terms in
    // conflicting ways (F4 on destruction vs Fr1 on production).
    // -----------------------------------------------------------------
    if (controlRCHellsten_ && controlRCMentor_)
    {
        FatalErrorInFunction
            << "Both controlRCHellsten (Hellsten 1998) and controlRCMentor "
            << "(Spalart-Shur / Menter-Smirnov 2009) are activated.\n"
            << "These curvature corrections are mutually exclusive.\n"
            << "Please set only ONE to 'yes' in the turbulenceProperties "
            << "dictionary."
            << exit(FatalError);
    }

    // Print active-extension summary at construction time
    if (stabilization_)
    {
        Info<< "    Larsen & Fuhrman (2018) stabilization: ACTIVE" << nl
            << "        lambda2  = " << lambda2_.value() << nl
            << "        alphaBS  = " << alphaBS_.value() << nl;
    }

    if (controlRCHellsten_)
    {
        Info<< "    Hellsten (1998) Rotation/Curvature Correction: ACTIVE" << nl
            << "        cRC      = " << cRC_.value() << nl;
    }

    if (controlRCMentor_)
    {
        Info<< "    Menter-Smirnov (2009) Curvature Correction: ACTIVE" << nl
            << "        cr1      = " << cr1_.value() << nl
            << "        cr2      = " << cr2_.value() << nl
            << "        cr3      = " << cr3_.value() << nl;
    }

    Info<< endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicEddyViscosityModel>
void kOmegaSSTNewBase<BasicEddyViscosityModel>::setDecayControl
(
    const dictionary& dict
)
{
    decayControl_.readIfPresent("decayControl", dict);

    if (decayControl_)
    {
        kInf_.read(dict);
        omegaInf_.read(dict);

        Info<< "    Employing decay control with kInf:" << kInf_
            << " and omegaInf:" << omegaInf_ << endl;
    }
    else
    {
        kInf_.value() = 0;
        omegaInf_.value() = 0;
    }
}


template<class BasicEddyViscosityModel>
bool kOmegaSSTNewBase<BasicEddyViscosityModel>::read()
{
    if (BasicEddyViscosityModel::read())
    {
        // Standard SST coefficients
        alphaK1_.readIfPresent(this->coeffDict());
        alphaK2_.readIfPresent(this->coeffDict());
        alphaOmega1_.readIfPresent(this->coeffDict());
        alphaOmega2_.readIfPresent(this->coeffDict());
        gamma1_.readIfPresent(this->coeffDict());
        gamma2_.readIfPresent(this->coeffDict());
        beta1_.readIfPresent(this->coeffDict());
        beta2_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());
        a1_.readIfPresent(this->coeffDict());
        b1_.readIfPresent(this->coeffDict());
        c1_.readIfPresent(this->coeffDict());
        F3_.readIfPresent("F3", this->coeffDict());

        setDecayControl(this->coeffDict());

        // Larsen & Fuhrman (2018) coefficients
        stabilization_.readIfPresent("stabilization", this->coeffDict());
        lambda2_.readIfPresent(this->coeffDict());
        alphaBS_.readIfPresent(this->coeffDict());

        // Hellsten (1998) RCH coefficients
        controlRCHellsten_.readIfPresent
        (
            "controlRCHellsten", this->coeffDict()
        );
        cRC_.readIfPresent(this->coeffDict());

        // Menter-Smirnov (2009) CC coefficients
        controlRCMentor_.readIfPresent
        (
            "controlRCMentor", this->coeffDict()
        );
        cr1_.readIfPresent(this->coeffDict());
        cr2_.readIfPresent(this->coeffDict());
        cr3_.readIfPresent(this->coeffDict());

        // Re-validate mutual exclusivity at run-time
        if (controlRCHellsten_ && controlRCMentor_)
        {
            FatalErrorInFunction
                << "Both controlRCHellsten and controlRCMentor are activated. "
                << "These are mutually exclusive."
                << exit(FatalError);
        }

        return true;
    }

    return false;
}


// =========================================================================
//  correct() — main time-step routine
//
//  Execution flow:
//
//  1.  Compute gradU → S2, pOmega, G, GbyNu0
//  2.  [STABILIZATION]  If active: pre-solve nut update with L&F limiter
//  3.  [BUOYANCY]       If rho + g found: compute tBuoyancy
//  4.  [HELLSTEN RCH]   If active: compute F4 from Richardson number
//  5.  [MENTER CC]      If active: compute Fr1 from DS_ij/Dt
//  6.  Assemble & solve omega equation  (F4 on destruction, Fr1 on prod.)
//  7.  Assemble & solve k equation      (Fr1 on production, buoyancy)
//  8.  End-of-step nut correction
//
//  When all extensions are OFF, the code path is identical to standard SST
//  from OpenFOAM v2512.
// =========================================================================
template<class BasicEddyViscosityModel>
void kOmegaSSTNewBase<BasicEddyViscosityModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // =================================================================
    // In practical terms, BasicEddyViscosityModel::correct() doesn't 
    // compute anything turbulence-specific — no field updates, 
    // no equation solving. 
    // It just ensures the base-class bookkeeping is done 
    // (flags set, counters updated) 
    // so that the infrastructure knows the turbulence model has been 
    // invoked this time step. 
    // All the real work (solving k, omega, updating nut) happens in 
    // the lines that follow in kOmegaSSTNewBase::correct().
    // =================================================================
    BasicEddyViscosityModel::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    // =================================================================
    // Step 1: Velocity gradient and invariants
    //
    // tgradU is computed ONCE.  From it we derive:
    //   S2      = 2*|S|²    (strain-rate invariant)
    //   pOmega  = 2*|W|²    (rotation-rate invariant)
    //
    // pOmega is needed by:
    //   - L&F stabilization limiter (when stabilization_ is ON)
    //   - Hellsten RCH correction   (when controlRCHellsten_ is ON)
    //   - Menter-Smirnov CC         (when controlRCMentor_ is ON)
    //
    // We allocate it only when at least one of these extensions needs
    // it.  On the purely-standard-SST path it is never computed.
    // =================================================================
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField S2(this->S2(tgradU()));

    // Determine whether pOmega is needed by any active extension
    const bool needPOmega =
        stabilization_ || controlRCHellsten_ || controlRCMentor_;

    tmp<volScalarField> tpOmega;
    if (needPOmega)
    {
        tpOmega = 2.0*magSqr(skew(tgradU()));
    }

    volScalarField::Internal GbyNu0(this->GbyNu0(tgradU(), S2));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);


    // =================================================================
    // Step 2: [STABILIZATION]  Pre-solve nut update
    //
    // When the L&F stabilization is active, nut is corrected BEFORE
    // the transport equations are assembled so that the implicit
    // diffusion terms (DkEff, DomegaEff) use the stabilized viscosity.
    // G is re-computed with the updated nut.
    // =================================================================
    if (stabilization_)
    {
        const volScalarField F23(this->F23());
        this->correctNut(S2, tpOmega(), F23);
        // Re-compute G with the updated (stabilized) nut
        G = volScalarField::Internal(this->GName(), nut*GbyNu0);
    }


    // =================================================================
    // Step 3: [BUOYANCY]  Dynamic detection of buoyancy production
    //
    // In multiphase VOF solvers (interFoam, interIsoFoam, etc.) two
    // objects live in the objectRegistry:
    //   - "rho"  (volScalarField)                  — mixture density
    //   - "g"    (uniformDimensionedVectorField)    — gravity
    //
    // When BOTH are present we form:
    //   Gb = alphaBS * (g · ∇ρ)          [kg/(m³·s²)]
    //
    // This is injected into the k-equation as:
    //   kEqn -= fvm::SuSp(alpha * nut * Gb / k, k)
    //
    // Sign analysis (z-up, water below air):
    //   g = (0,0,-9.81), grad(rho)_z < 0
    //   → g·∇ρ > 0  (stable)  → implicit SINK  → k suppressed  ✓
    //   Unstable: g·∇ρ < 0    → explicit SOURCE → k enhanced    ✓
    //
    // In single-phase runs tBuoyancy stays invalid → no overhead.
    // =================================================================
    /*
    tmp<volScalarField> tBdensity;

    if (this->db().template foundObject<volScalarField>("rho") &&
        this->db().template foundObject<uniformDimensionedVectorField>("g"))
    {
        const volScalarField& varRho =
            this->db().template
            lookupObject<volScalarField>("rho");

        const uniformDimensionedVectorField& g =
            this->db().template
            lookupObject<uniformDimensionedVectorField>("g");

        // tBdensity = alphaBS * (g · ∇ρ)
        // Dimensions: [m/s^2] × [kg/m^4] = [kg/(m^3·s^2)] = [N/m^5]
        // nut·tBdensity/k → [1/s]: correct dimensions for fvm::SuSp.
        tBdensity = alphaBS_ * this->nut_ * (g & fvc::grad(varRho));
    }
    */
    tmp<volScalarField> tBdensity;

    if
    (
        stabilization_
     && this->mesh_.template
            foundObject<volScalarField>("rho")
     && this->mesh_.template
            foundObject<uniformDimensionedVectorField>("g")
    )
    {
        const volScalarField& varRho =
            this->mesh_.template
                lookupObject<volScalarField>("rho");

        const uniformDimensionedVectorField& g =
            this->mesh_.template
                lookupObject<uniformDimensionedVectorField>("g");

        tBdensity = alphaBS_ * (g & fvc::grad(varRho));
    }


    // =================================================================
    // Step 4: [HELLSTEN RCH]  Rotation/Curvature Correction
    //
    //   F4 = 1 / (1 + cRC * Ri)
    //
    //   where the Richardson number is:
    //     Ri = (W/S) * (W/S - 1)
    //     W  = sqrt(pOmega) = sqrt(2*Ω_ij*Ω_ij)   — rotation rate
    //     S  = sqrt(S2)     = sqrt(2*S_ij*S_ij)     — strain rate
    //
    //   Physical interpretation:
    //     - Pure shear (W = S):       Ri = 0  → F4 = 1     (no effect)
    //     - Stabilising rotation W>S: Ri > 0  → F4 < 1
    //       → less ω destruction → more ω → less k  (suppresses turbulence)
    //     - Destabilising curvature W<S: Ri < 0  → F4 > 1
    //       → more ω destruction → less ω → more k (enhances turbulence)
    //
    //   F4 multiplies the omega destruction term:  -Sp(F4*beta*omega, omega)
    //   When RCH is OFF, tF4 = 1.0 everywhere (no modification).
    // =================================================================
    tmp<volScalarField::Internal> tF4
    (
        tmp<volScalarField::Internal>::New
        (
            IOobject
            (
                IOobject::scopedName(this->type(), "F4"),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("1", dimless, 1.0)
        )
    );

    if (controlRCHellsten_)
    {
        // Small value to prevent division by zero where S → 0
        dimensionedScalar vSmall("VSMALL", dimless/dimTime, VSMALL);

        // Compute rotation rate W = sqrt(2*Ω_ij*Ω_ij) and
        //         strain rate  S = sqrt(2*S_ij*S_ij)
        // from the pre-computed invariants pOmega and S2.
        const volScalarField::Internal W(sqrt(tpOmega()()));
        const volScalarField::Internal S(sqrt(S2()));

        // Richardson number: Ri = (W/S) * (W/S - 1)
        //   Ri = 0  when W = S  (simple shear)
        //   Ri > 0  when W > S  (stabilising rotation)
        //   Ri < 0  when W < S  (destabilising curvature)
        const volScalarField::Internal Ri
        (
            W/max(S, vSmall) * (W/max(S, vSmall) - 1.0)
        );

        // F4 = 1/(1 + cRC * Ri)
        tF4.ref() = 1.0 / (1.0 + cRC_.value() * Ri);
    }


    // =================================================================
    // Step 5: [MENTER CC]  Spalart-Shur / Menter-Smirnov Curvature
    //                       Correction
    //
    //   Fr1 = max( min( f_rotation, 1.25 ), 0.0 )
    //
    //   where:
    //     f_rotation = (1 + cr1) * 2r*/(1 + r*)
    //                  * (1 - cr3*atan(cr2*r̂)) - cr1
    //
    //     r*    = |S|/|W| = sqrt(S2/pOmega)  (strain-to-rotation ratio)
    //     D     = sqrt(max(S2, 0.09*ω²))     (combined scale, ensures D>0
    //                                          even in irrotational regions)
    //     r̂     = 2*(W·S)::(DS/Dt) / (|W|*D³) (Lagrangian curvature measure)
    //
    //   DS_ij/Dt is the material derivative of the strain-rate tensor,
    //   computed as ddt(S) + div(alphaRhoPhi, S).  This is what makes
    //   the correction "Lagrangian" — it tracks how the strain field
    //   evolves along streamlines, which is the physical mechanism by
    //   which streamline curvature affects turbulence.
    //
    //   Fr1 multiplies the production terms in both equations:
    //     Pω  → gamma*GbyNu0*Fr1
    //     Pk  → Pk*Fr1
    //   When CC is OFF, tFr1 = 1.0 everywhere (no modification).
    //
    //   Bounds: Fr1 ∈ [0, 1.25]
    //     Fr1 > 1  → enhanced production (destabilising curvature)
    //     Fr1 < 1  → suppressed production (stabilising rotation)
    //     Fr1 = 1  → no effect (simple shear or irrotational flow)
    // =================================================================
    tmp<volScalarField::Internal> tFr1
    (
        tmp<volScalarField::Internal>::New
        (
            IOobject
            (
                IOobject::scopedName(this->type(), "Fr1"),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("1", dimless, 1.0)
        )
    );

    if (controlRCMentor_)
    {
        // Symmetric and skew-symmetric parts of gradU (already available
        // from tgradU which is still alive at this point)
        const volSymmTensorField S(symm(tgradU()));
        const volTensorField     W(skew(tgradU()));

        // r* = |S|/|W| = sqrt(S2/pOmega)
        // Clamped from below to avoid 0/0 where W → 0.
        const volScalarField::Internal rStar
        (
            sqrt
            (
                S2()
              / max
                (
                    tpOmega()(),
                    dimensionedScalar
                    (
                        "minW2", 
                        tpOmega().dimensions(), 
                        SMALL
                    )
                )
            )
        );

        // D = sqrt(max(S2, 0.09*omega²))
        // The 0.09*omega² floor ensures D remains positive in
        // irrotational regions where S → 0.
        const volScalarField::Internal D
        (
            sqrt
            (
                max
                (
                    S2(), 
                    0.09*sqr(omega_())
                )
            )
        );

        // |W| = sqrt(pOmega), clamped from below
        const volScalarField::Internal Wmag
        (
            sqrt
            (
                max
                (
                    tpOmega()(),
                    dimensionedScalar
                    (
                        "minW2", 
                        tpOmega().dimensions(), 
                        SMALL
                    )
                )
            )
        );

        // DS_ij/Dt: Lagrangian (material) derivative of the strain tensor.
        // Computed as ddt(S) + div(alphaRhoPhi, S).
        // Note: this uses the face mass flux for the convective part.
        const volSymmTensorField DSiDt
        (
            fvc::ddt(S) + fvc::div(alphaRhoPhi, S)
        );

        // r̂ = 2*(W·S)::(DS/Dt) / (|W| * D³)
        //
        // (W & S)_ij = W_ik * S_kj = W_ik * S_jk  (S is symmetric)
        // The double contraction && with DS_ij/Dt gives the scalar
        // numerator of the curvature measure.
        const volScalarField::Internal rTilda
        (
            2.0*((W() & S()) && DSiDt()) / (Wmag*pow3(D))
        );

        // f_rotation = (1+cr1) * 2r*/(1+r*) * (1 - cr3*atan(cr2*r̂)) - cr1
        // Fr1 = max(min(f_rotation, 1.25), 0.0)
        //
        // The atan() term saturates r̂ so that extreme curvature
        // does not cause unbounded correction.
        tFr1.ref() = max
        (
            min
            (
                (1.0 + cr1_.value())*2.0*rStar/(1.0 + rStar)
               *(1.0 - cr3_.value()*atan(cr2_.value()*rTilda))
               - cr1_.value(),
                1.25
            ),
            0.0
        );
    }


    // =================================================================
    // Step 6: Omega boundary condition update + equation assembly
    // =================================================================

    // - boundary condition changes a cell value
    // - normally this would be triggered through correctBoundaryConditions
    // - which would do
    //      - fvPatchField::evaluate() which calls
    //      - fvPatchField::updateCoeffs()
    // - however any processor boundary conditions already start sending
    //   at initEvaluate so would send over the old value.
    // - avoid this by explicitly calling updateCoeffs early and then
    //   only doing the boundary conditions that rely on initEvaluate
    //   (currently only coupled ones)

    //- 1. Explicitly swap values on coupled boundary conditions
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();
    // omegaWallFunctions change the cell value! Make sure to push these to
    // coupled neighbours. Note that we want to avoid the re-updateCoeffs
    // of the wallFunctions so make sure to bypass the evaluate on
    // those patches and only do the coupled ones.
    omega_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();
    ////- 2. Make sure the boundary condition calls updateCoeffs from
    ////     initEvaluate
    ////     (so before any swap is done - requires all coupled bcs to be
    ////      after wall bcs. Unfortunately this conflicts with cyclicACMI)
    //omega_.correctBoundaryConditions();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));
    const volScalarField F23(this->F23());

    {
        const volScalarField::Internal gamma(this->gamma(F1));
        const volScalarField::Internal beta(this->beta(F1));

        // Apply the Bradshaw / stagnation-point limiter to GbyNu0
        GbyNu0 = GbyNu(GbyNu0, F23(), S2());

        // --- Omega equation ---
        //
        // Key modifications vs standard SST:
        //   Production:  gamma*GbyNu0 * tFr1
        //     (tFr1 = 1.0 unless Menter CC is active)
        //   Destruction: tF4 * beta*omega²
        //     (tF4 = 1.0 unless Hellsten RCH is active)
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
         ==
            // Production: gamma * (G/nu) * Fr1
            alpha()*rho()*gamma*GbyNu0*tFr1()
            // Compressibility correction
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
            // Destruction: F4 * beta * omega²
            // F4 modulates destruction under rotation/curvature effects
          - fvm::Sp(tF4()*alpha()*rho()*beta*omega_(), omega_)
            // SST cross-diffusion: active in outer (k-ε) layer only
            // (F1-1) = 0 near wall, -1 in free stream
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
            // Far-field decay source (zero if decayControl=no)
          + alpha()*rho()*beta*sqr(omegaInf_)
          + Qsas(S2(), gamma, beta)
          + omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    {
        // --- k equation ---
        //
        // Key modifications vs standard SST:
        //   Production:  Pk(G) * tFr1
        //     (tFr1 = 1.0 unless Menter CC is active)
        //   Buoyancy:    -= SuSp(alpha*nut*Gb/k, k)
        //     (only when tBdensity is valid, i.e. multiphase)
        tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(alpha, rho, k_)
          + fvm::div(alphaRhoPhi, k_)
          - fvm::laplacian(alpha*rho*DkEff(F1), k_)
         ==
            // Production * Fr1
            alpha()*rho()*Pk(G)*tFr1()
            // Compressibility correction
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
            // Destruction: betaStar * omega * k
          - fvm::Sp(alpha()*rho()*epsilonByk(F1, tgradU()), k_)
            // Far-field decay source (zero if decayControl=no)
          + alpha()*rho()*betaStar_*omegaInf_*kInf_
          + kSource()
          + fvOptions(alpha, rho, k_)
        );

        // [BUOYANCY]  Inject buoyancy production into k-equation.
        //
        // "-= SuSp(coeff/k, k)" means:
        //   coeff > 0  → implicit sink   (stable stratification)
        //   coeff < 0  → explicit source  (unstable stratification)
        //
        // k denominator is guarded with kMin_ to prevent division by
        // zero in regions where k → 0 (e.g. far-field air).
        if (tBdensity.valid())
        {
            kEqn.ref() -= fvm::SuSp
            (
                // alpha() * 
                this->nut_() * tBdensity()() / max(k_(), this->kMin_),
                k_
            );
        }

        // tgradU is no longer needed; release it before solve.
        tgradU.clear();

        kEqn.ref().relax();
        fvOptions.constrain(kEqn.ref());
        solve(kEqn);
        fvOptions.correct(k_);
        bound(k_, this->kMin_);
    }

    // =================================================================
    // Step 8: End-of-step nut correction
    //
    // Routes to the stabilized 3-arg path or the standard 1-arg path
    // depending on the stabilization switch.
    // =================================================================
    if (stabilization_)
    {
        correctNut(S2, tpOmega(), F23);
    }
    else
    {
        correctNut(S2);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

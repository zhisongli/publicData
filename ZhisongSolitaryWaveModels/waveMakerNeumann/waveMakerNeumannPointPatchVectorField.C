/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 IH-Cantabria
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

#include "waveMakerNeumannPointPatchVectorField.H"
#include "mathematicalConstants.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "gravityMeshObject.H"

#include "polyMesh.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum<Foam::waveMakerNeumannPointPatchVectorField::motionTypes>
Foam::waveMakerNeumannPointPatchVectorField::motionTypeNames
({
    { motionTypes::piston, "piston" },
    { motionTypes::flap, "flap" },
    { motionTypes::solitary, "solitary" },
    { motionTypes::solitaryBoussinesq, "solitaryBoussinesq" },
    { motionTypes::solitaryWeakly, "solitaryWeakly" },
    { motionTypes::solitaryStrongly, "solitaryStrongly" },
    { motionTypes::solitaryCF2003, "solitaryCF2003" }
});




// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

const Foam::vector& Foam::waveMakerNeumannPointPatchVectorField::g()
{
    const meshObjects::gravity& gf = meshObjects::gravity::New(db().time());

    if (mag(gf.value()) < SMALL)
    {
        FatalErrorInFunction
            << "Gravity vector is not set.  Please update "
            << gf.uniformDimensionedVectorField::path()
            << exit(FatalError);
    }

    return gf.value();
}


Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::waveLength
(
    const scalar h,
    const scalar T
)
{
    const scalar L0 = mag(g())*T*T/(constant::mathematical::twoPi);
    scalar L = L0;

    for (label i=1; i<=100; ++i)
    {
        L = L0*tanh(constant::mathematical::twoPi*h/L);
    }

    return L;
}


Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::timeCoeff
(
    const scalar t
) const
{
    return clamp(t/rampTime_, zero_one{});
}


void Foam::waveMakerNeumannPointPatchVectorField::initialiseGeometry()
{
    // Global patch extents
    const vectorField& Cp = this->patch().localPoints();
    const vectorField CpLocal(Cp);
    boundBox bb(CpLocal, true);

    const scalar xMin = bb.min().x();
    const scalar xMax = bb.max().x();
    const scalar yMin = bb.min().y();
    const scalar yMax = bb.max().y();
    zSpan_ = bb.max().z() - bb.min().z();

    zMinGb_ = bb.min().z();
    reduce(zMinGb_, minOp<scalar>());

    // Global x, y positions of the paddle centres
    xPaddle_.setSize(nPaddle_, 0);
    yPaddle_.setSize(nPaddle_, 0);
    const scalar xMid = xMin + 0.5*(xMax - xMin);
    const scalar paddleDy = (yMax - yMin)/scalar(nPaddle_);

    for (label paddlei = 0; paddlei < nPaddle_; ++paddlei)
    {
        xPaddle_[paddlei] = xMid;
        yPaddle_[paddlei] = paddlei*paddleDy + yMin + 0.5*paddleDy;
    }

    // Local face centres
    x_ = this->patch().localPoints().component(0);
    y_ = this->patch().localPoints().component(1);
    z_ = this->patch().localPoints().component(2);

    // Local point-to-paddle addressing
    pointToPaddle_.setSize(this->patch().size(), -1);

    forAll(pointToPaddle_, ppi)
    {
        pointToPaddle_[ppi] = floor((y_[ppi] - yMin)/(paddleDy+0.01*paddleDy));
    }
}

#include "_funGammaEpsilon.C"
#include "_funKappa.C"
#include "_funCelerity.C"
#include "_funPeriod.C"
#include "_funStroke.C"
#include "_funEta.C"
#include "_funDXiDt.C"

#include "_funCF2003.C"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveMakerNeumannPointPatchVectorField::waveMakerNeumannPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF),
    motionType_(motionTypes::piston),
    n_(Zero),
    gHat_(Zero),
    initialDepth_(0),
    wavePeriod_(0),
    waveHeight_(0),
    wavePhase_(0),
    waveAngle_(0),
    startTime_(0),
    rampTime_(1),
    order_(0),
    Goring0MMT1_(0),
    secondOrder_(false),
    nPaddle_(0)
{}


Foam::waveMakerNeumannPointPatchVectorField::waveMakerNeumannPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict, IOobjectOption::NO_READ),
    motionType_(motionTypeNames.get("motionType", dict)),
    n_(dict.get<vector>("n")),
    gHat_(Zero),
    initialDepth_(dict.get<scalar>("initialDepth")),
    wavePeriod_(dict.get<scalar>("wavePeriod")),
    waveHeight_(dict.get<scalar>("waveHeight")),
    wavePhase_(dict.get<scalar>("wavePhase")),
    waveAngle_(dict.getOrDefault<scalar>("waveAngle", 0)),
    startTime_
    (
        dict.getOrDefault<scalar>
        (
            "startTime",
            db().time().startTime().value()
        )
    ),
    rampTime_(dict.get<scalar>("rampTime")),
    order_(dict.get<scalar>("order")),
    Goring0MMT1_(dict.get<scalar>("Goring0MMT1")),
    secondOrder_(dict.getOrDefault<bool>("secondOrder", false)),
    nPaddle_(dict.getOrDefault<label>("nPaddle", 1))
{
    // Create the co-ordinate system
    if (mag(n_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Patch normal direction vector is not set. 'n' = " << n_
            << exit(FatalIOError);
    }
    n_.normalise();

    gHat_ = (g() - n_*(n_&g()));
    if (mag(gHat_) < SMALL)
    {
        FatalIOErrorInFunction(dict)
            << "Patch normal and gravity directions must not be aligned. "
            << "'n' = " << n_ << " 'g' = " << g()
            << exit(FatalIOError);
    }
    gHat_.normalise();

    waveAngle_ *= constant::mathematical::pi/180;

    initialiseGeometry();

    waterDepthRef_.setSize(nPaddle_, -1);

    if (!dict.found("value"))
    {
        updateCoeffs();
    }

    k_nu_A_kk_C_nuk_u_uw_UnitDepth = k_nu_A_kk_C_nuk_u_uw(waveHeight_ / initialDepth_, initialDepth_ / initialDepth_);
    
}


Foam::waveMakerNeumannPointPatchVectorField::waveMakerNeumannPointPatchVectorField
(
    const waveMakerNeumannPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveAngle_(ptf.waveAngle_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    order_(ptf.order_),
    Goring0MMT1_(ptf.Goring0MMT1_),
    secondOrder_(ptf.secondOrder_),
    nPaddle_(ptf.nPaddle_)
{}


Foam::waveMakerNeumannPointPatchVectorField::waveMakerNeumannPointPatchVectorField
(
    const waveMakerNeumannPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    motionType_(ptf.motionType_),
    n_(ptf.n_),
    gHat_(ptf.gHat_),
    initialDepth_(ptf.initialDepth_),
    wavePeriod_(ptf.wavePeriod_),
    waveHeight_(ptf.waveHeight_),
    wavePhase_(ptf.wavePhase_),
    waveAngle_(ptf.waveAngle_),
    startTime_(ptf.startTime_),
    rampTime_(ptf.rampTime_),
    order_(ptf.order_),
    Goring0MMT1_(ptf.Goring0MMT1_),
    secondOrder_(ptf.secondOrder_),
    nPaddle_(ptf.nPaddle_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveMakerNeumannPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (firstTime == 0)
    {
        // Set the reference water depth
        if (initialDepth_ != 0 )
        {
            forAll(waterDepthRef_, paddlei)
            {
                waterDepthRef_[paddlei] = initialDepth_;
            }
        }
        else
        {
            FatalErrorInFunction
               << "initialDepth is not set.  Please update "
               << abort(FatalError);
        }


        Info<< " WaterDepth at the wavepaddles = " << waterDepthRef_ << endl;
        firstTime = 1;
    }

    const scalar t = db().time().value() - startTime_;

    scalarField waveLength_(nPaddle_, -1);

    scalarField waveK(nPaddle_, -1);
    scalarField waveKx(nPaddle_, -1);
    scalarField waveKy(nPaddle_, -1);

    forAll(waveK, padddlei)
    {
        waveLength_[padddlei] =
            waveLength(waterDepthRef_[padddlei], wavePeriod_);

        waveK[padddlei] = constant::mathematical::twoPi/waveLength_[padddlei];
        waveKx[padddlei] = waveK[padddlei]*cos(waveAngle_);
        waveKy[padddlei] = waveK[padddlei]*sin(waveAngle_);
    }
    const scalar sigma = 2*constant::mathematical::pi/wavePeriod_;

    switch (motionType_)
    {

        #include "srcFlapDefault.C"

        #include "srcPistonDefault.C"
        
        #include "srcSolitaryDefault.C"

        #include "srcBoussinesq.C"

        #include "srcStrongly.C"
        #include "srcWeakly.C"

        #include "srcCF2003.C"
        
        default:
        {
            FatalErrorInFunction
                << "Unhandled enumeration " << motionTypeNames[motionType_]
                << abort(FatalError);
        }
    }

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void Foam::waveMakerNeumannPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    os.writeEntry("motionType", motionTypeNames[motionType_]);
    os.writeEntry("n", n_);
    os.writeEntry("initialDepth", initialDepth_);
    os.writeEntry("wavePeriod", wavePeriod_);
    os.writeEntry("waveHeight", waveHeight_);
    os.writeEntry("wavePhase", wavePhase_);
    os.writeEntry("waveAngle", waveAngle_);
    os.writeEntry("startTime", startTime_);
    os.writeEntry("rampTime", rampTime_);
    os.writeEntry("order", order_);
    os.writeEntry("Goring0MMT1", Goring0MMT1_);
    os.writeEntry("secondOrder", secondOrder_);
    os.writeEntry("nPaddle", nPaddle_);

    this->writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePointPatchTypeField
    (
        pointPatchVectorField,
        waveMakerNeumannPointPatchVectorField
    );
}

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 IH-Cantabria
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "Clamond1999WaveModel.H"
#include "addToRunTimeSelectionTable.H"
// #include "complex.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(Clamond1999, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        Clamond1999,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::Clamond1999::waveCelerity
(
    const scalar H,
    const scalar h
) const
{
    return sqrt(mag(g_)*(H + h));
}

Foam::scalar Foam::waveModels::Clamond1999::halfLength
(
    const scalar H,
    const scalar h
) const
{
    return 3.5*h/sqrt(H/h);
}

// sech
Foam::scalar Foam::waveModels::Clamond1999::sech
(
    const scalar x
) const
{
    return 1.0 / cosh(x);
}

Foam::complex Foam::waveModels::Clamond1999::sech
(
    const complex x
) const
{
    return 1.0 / cosh(x);
}

// csc
Foam::scalar Foam::waveModels::Clamond1999::csc
(
    const scalar x
) const
{
    scalar small = SMALL;
    return 1.0 / (sin(x) + small);
}

Foam::complex Foam::waveModels::Clamond1999::csc
(
    const complex x
) const
{
    scalar small = SMALL;
    return 1.0 / (sin(x) + small);
}

// sec
Foam::scalar Foam::waveModels::Clamond1999::sec
(
    const scalar x
) const
{
    scalar small = SMALL;
    return 1.0 / (cos(x) + small);
}

Foam::complex Foam::waveModels::Clamond1999::sec
(
    const complex x
) const
{
    scalar small = SMALL;
    return 1.00 / (cos(x) + small);
}


Foam::vector Foam::waveModels::Clamond1999::kAC
(
    const scalar H,
    const scalar h
) const
{
    const scalar g = mag(g_);
    scalar epsilon = H / h;
    scalar a = H;

    complex inK = complex(1.0 / (2.0 * h), 0.0);
    complex inC = complex(sqrt(g * h) * 0.9, 0.0);
    complex inA = 2.0 * epsilon * inC / 2.0;

    // tmpY
    complex sZero(complex(Zero, Zero));
    Vector<complex> tmpY(sZero + inK, sZero + inA, sZero + inC);

    label N(10);
    label iter(0);
    const scalar eps(1.e-6);
    scalar residual(1.0);
    scalar residualB(1.0);
    scalar residualC(1.0);


    Info << "Initial guess of tmpY = " << tmpY << endl;

    while (iter <= N)
    {
        Info << "===================================================" << endl;
        Info << "iter = " << iter << endl;

        Info << "Initial tmpY = " << tmpY << endl;

        iter++;

        complex k = complex(inK);
        complex A = complex(inA);
        complex C = complex(inC);

        complex kh = k * h;

        complex C2kh = csc(2*kh);
        complex S1kh = sec(kh);
        complex S2kh = sec(2*kh);
        complex T1kh = tan(kh);
        complex T2kh = tan(2*kh);
        complex S1kah = sec(k * (a + h));
        complex T1kah = tan(k*(a + h));

        // composing the Jacobian matrix
        complex JacobianM11 = a*C - A*pow(S1kh,2)*(a + h + 2*a*kh*T1kh);
        complex JacobianM12 = -(a*k*pow(S1kh,2)) - T1kh;
        complex JacobianM13 = a*k;

        complex JacobianM21 = -2*A*(a + h)*pow(S1kah,2)*T1kah;
        complex JacobianM22 = -pow(S1kah,2);
        complex JacobianM23 = 1 - C / sqrt(complex(pow(C,2) - 2*a*g));

        complex JacobianM31 =(C*h)/(sqrt(2.0)* sqrt(kh)) - C2kh*h*S2kh*sqrt(g*h*T2kh);
        complex JacobianM32 = complex(0.0);
        complex JacobianM33 = sqrt(2.0) * sqrt(kh);

        Tensor<complex> JacobianM = {
            {JacobianM11, JacobianM12, JacobianM13},
            {JacobianM21, JacobianM22, JacobianM23},
            {JacobianM31, JacobianM32, JacobianM33}
        };

        // Info << "JacobianM = " << JacobianM << endl;

        //
        complex rightB1 = a*k*(C - A*pow(S1kh,2)) - A*T1kh;
        complex rightB2 = C - sqrt(complex(pow(C,2) - 2*a*g)) - A*pow(S1kah,2);
        complex rightB3 = sqrt(2.0)*C*sqrt(kh) - sqrt(g*h*T2kh);
        Vector<complex> rightB = {rightB1, rightB2, rightB3};

        // Info << "rightB = " << rightB << endl;

        tmpY -= (inv(JacobianM) & rightB);

        Info << "updated tmpY = " << tmpY << endl;


        residual = Foam::mag(inK - tmpY[0]);
        residualB = Foam::mag(inA - tmpY[1]);
        residualC = Foam::mag(inC - tmpY[2]);
        inK = tmpY[0];
        inA = tmpY[1];
        inC = tmpY[2];

        // Info << "inK = " << inK << endl;
        // Info << "residual = " << residual << endl;
        
        if ((residual < eps) && (residualB < eps) && (residualC < eps))
        {
            Info << "===================================================" << endl;
            Info << "kAC are coverged after " << iter << " iteration(s)." << endl;
            Info << "Wave height = " << H << " m" << endl;
            Info << "Water depth = " << h << " m" << endl;
            Info << "      kappa = " << inK << endl;
            Info << "          A = " << inA << endl;
            Info << "          C = " << inC << endl;
            Info << "===================================================" << endl;
            return vector(inK.real(), inA.real(), inC.real());
        }
        else if (iter > N-1)
        {
            FatalErrorInFunction
                << "Newton-Raphson iterations diverging: "
                << "iterations = " << iter
                << exit(FatalError);
        }
       
    }
    
    // return vector(0.0, 0.0, 0.0);
}


Foam::scalar Foam::waveModels::Clamond1999::eta
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0
) const
{
    // vector kAC = this->kAC(H, h);
    scalar k = kAC_[0];
    scalar A = kAC_[1];
    scalar C = kAC_[2];

    // scalar ts = 3.5*h/sqrt(H/h);
    scalar ts = this->halfLength(H, h);
    scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    scalar Xb = Xa / 1.0;

    scalar kh = (k * h);
    
    scalar tmpA1 = A*sin(2*kh);
    scalar tmpA2 = C*k*(cos(2*kh) + cosh(2*k*Xb));
    scalar tmpA = tmpA1 / tmpA2;

    scalar tmpB1 = 2.0*A*(1.0 + cos(2*kh)*cosh(2*k*Xb));
    scalar tmpB2 = C*pow(cos(2*kh) + cosh(2*k*Xb),2);
    scalar tmpB = 1.0 - tmpB1 / tmpB2;
    
    return tmpA / tmpB;
}




Foam::vector Foam::waveModels::Clamond1999::Uf
(
    const scalar H,
    const scalar h,
    const scalar x,
    const scalar y,
    const scalar theta,
    const scalar t,
    const scalar X0,
    const scalar z
) const
{
    // vector kAC = this->kAC(H, h);
    scalar k = kAC_[0];
    scalar A = kAC_[1];
    scalar C = kAC_[2];


    // const scalar ts = 3.5*h/sqrt(H/h);
    const scalar ts = this->halfLength(H, h);
    const scalar Xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    const scalar Xb = Xa / 1.0;

    // const scalar kh = (k * h);
    const scalar Sech1kX = sech(k * Xb);
    const scalar Tanh1kX = tanh(k * Xb);
    const scalar Sin1kZ = sin(k * z);
    const scalar Sin2kZ = sin(2 * k * z);
    const scalar Cos2kZ = cos(2*k*z);

    const scalar w = (A*pow(Sech1kX,2)*Sin2kZ*Tanh1kX)/pow(1 - pow(Sech1kX,2)*pow(Sin1kZ,2),2);
    
    scalar u = (A*(Cos2kZ*pow(Sech1kX,2) + pow(Sech1kX,4)*pow(Sin1kZ,2)))/pow(1 - pow(Sech1kX,2)*pow(Sin1kZ,2),2);
    
    const scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    return vector(u, v, w);
}


void Foam::waveModels::Clamond1999::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    forAll(level, paddlei)
    {
        const scalar eta =
            this->eta
            (
                waveHeight_,
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveAngle_,
                t,
                x0_
            );

        level[paddlei] = waterDepthRef_ + tCoeff*eta;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::Clamond1999::Clamond1999
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    solitaryWaveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        // Info << "Clamond1999 reading file..." << endl;
        readDict(dict);
    }
    Info << "   kAC will be calculated." << endl;
    // Info << "   waveHeight_ = " << waveHeight_ << endl;
    // Info << "   waterDepthRef_ = " << waterDepthRef_ << endl;

    kAC_ = this->kAC(waveHeight_, waterDepthRef_);
    // Info << "   kAC_k = " << kAC_[0] << endl;
    // Info << "   kAC_A = " << kAC_[1] << endl;
    // Info << "   kAC_C = " << kAC_[2] << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::Clamond1999::readDict(const dictionary& overrideDict)
{
    if (solitaryWaveModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::Clamond1999::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    forAll(U_, facei)
    {
        // Fraction of geometry represented by paddle - to be set
        scalar fraction = 1;

        // Height - to be set
        scalar z = 0;

        setPaddlePropeties(level, facei, fraction, z);

        if (fraction > 0)
        {
            const label paddlei = faceToPaddle_[facei];

            const vector Uf = this->Uf
            (
                waveHeight_,
                waterDepthRef_,
                xPaddle_[paddlei],
                yPaddle_[paddlei],
                waveAngle_,
                t,
                x0_,
                z
            );

            U_[facei] = fraction*Uf*tCoeff;
        }
    }
}


void Foam::waveModels::Clamond1999::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //

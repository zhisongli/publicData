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

#include "GrimshawNewWaveModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(GrimshawNew, 0);
    addToRunTimeSelectionTable
    (
        waveModel,
        GrimshawNew,
        patch
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::scalar Foam::waveModels::GrimshawNew::alfa
(
    const scalar H,
    const scalar h
) const
{
    scalar eps = H/h;
    return sqrt(0.75*eps)*(1.0 - 0.625*eps + (71.0/128.0)*eps*eps);
    // scalar epsilon = H/h;
    // return sqrt(3.0*epsilon/4.0)*( 1.0 - 5.0/8.0*epsilon + 71.0/128.0*pow(epsilon,2) );
}

Foam::scalar Foam::waveModels::GrimshawNew::celerity
(
    const scalar H,
    const scalar h
) const
{
    const scalar G = mag(g_);
    scalar epsilon = H/h;
    return sqrt(G*h)*sqrt(1.0 + epsilon - 1.0/20.0*pow(epsilon,2) - 3.0/70.0*pow(epsilon,3));
}

Foam::scalar Foam::waveModels::GrimshawNew::waveLength
(
    const scalar H,
    const scalar h
) const
{
    return 2.0*4.22*h/sqrt(H/h);
}

Foam::scalar Foam::waveModels::GrimshawNew::wavePeriod
(
    const scalar H,
    const scalar h
) const
{
    return this->waveLength(H, h) / this->celerity(H, h);
}


Foam::scalar Foam::waveModels::GrimshawNew::eta
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
    // // olaFlow
    // double epsilon = H/h;
    // double beta = sqrt(3.0*epsilon/4.0)*( 1.0 - 5.0/8.0*epsilon
    //     + 71.0/128.0*pow(epsilon,2) );
    // double C = this->celerity(H, h);

    // double ts = this->waveLength(H, h)/2.0;
    // double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

    // double s = 1.0/cosh(beta/h*Xa);
    // double q = tanh(beta/h*Xa);

    // double sup = pow(s,2) - 3.0/4.0*epsilon*pow(s,2)*pow(q,2)
    //     + pow(epsilon,2)*( 5.0/8.0*pow(s,2)*pow(q,2) - 101.0/80.0*pow(s,4)*pow(q,2) );

    // return H*sup;

    // original
    const scalar eps = H/h;
    const scalar eps2 = eps*eps;
    const scalar eps3 = eps*eps2;

    // const scalar C =
    //     sqrt(mag(g_)*h)*sqrt(1.0 + eps - 0.05*eps2 - (3.0/70.0)*eps3);
    const scalar C = this->celerity(H, h);

    // const scalar ts = 3.5*h/sqrt(H/h);
    const scalar ts = this->waveLength(H, h) / 2.0;
    const scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    const scalar alfa = this->alfa(H, h);

    const scalar s = (1.0)/(cosh(alfa*(xa/h)));
    const scalar s2 = s*s;
    const scalar q = tanh(alfa*(xa/h));
    const scalar q2 = q*q;

    return h*(
            eps*s2
              - 0.75*eps2*s2*q2
              + eps3*(0.625*s2*q2 - 1.2625*s2*s2*q2)
            );
    
}


Foam::vector Foam::waveModels::GrimshawNew::Uf
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
    // const scalar G = mag(g_);

    // double epsilon = H/h;
    // double beta = sqrt(3.0*epsilon/4.0)*( 1.0 - 5.0/8.0*epsilon
    //     + 71.0/128.0*pow(epsilon,2) );
    // double C = this->celerity(H, h);

    // double ts = this->waveLength(H, h)/2.0;
    // double Xa = -C * t + ts - X0 + x * cos(theta) + y * sin(theta);

    // double s = 1.0/cosh(beta/h*Xa);
    // double q = tanh(beta/h*Xa);

    // // u, v
    // double u = pow(s,2)*epsilon
    //     - pow(epsilon,2)*
    //     (
    //         -1.0/4.0*pow(s,2) 
    //         + pow(s,4) 
    //         + pow(z/h,2)*(3.0/2.0*pow(s,2)-9.0/4.0*pow(s,4)) 
    //     )
    //     - pow(epsilon,3)*
    //     ( 
    //         19.0/40.0*pow(s,2) + 1.0/5.0*pow(s,4) - 6.0/5.0*pow(s,6)
    //         + pow(z/h,2)*( -3.0/2.0*pow(s,2)-15.0/4.0*pow(s,4)+15.0/2.0*pow(s,6) )
    //         + pow(z/h,4)*( -3.0/8.0*pow(s,2)+45.0/16.0*pow(s,4)-45.0/16.0*pow(s,6) ) 
    //     );

    // scalar v = u*sin(waveAngle_);
    // u *= cos(waveAngle_);

    // // w
    // double w = -sqrt(3.0*epsilon)*(z/h)*q*
    // (
    //     - epsilon*pow(s,2)
    //     + pow(epsilon,2)*
    //     ( 
    //         3.0/8.0*pow(s,2) 
    //         + 2.0*pow(s,4) 
    //         + pow(z/h,2)*(1.0/2.0*pow(s,2) - 3.0/2.0*pow(s,4) ) 
    //     )
    //     + pow(epsilon,3)*
    //     ( 
    //         49.0/640.0*pow(s,2) - 17.0/20.0*pow(s,4) - 18.0/5.0*pow(s,6) 
    //         + pow(z/h,2)*(-13.0/16.0*pow(s,2) - 25.0/16.0*pow(s,4) + 15.0/2.0*pow(s,6) )
    //         + pow(z/h,4)*(-3.0/40.0*pow(s,2)+9.0/8.0*pow(s,4)-27.0/16.0*pow(s,6) ) 
    //     ) 
    // );
    

    // u *= sqrt(G*h);
    // v *= sqrt(G*h);
    // w *= sqrt(G*h);

    // return vector(u, v, w);
    
    
    const scalar eps = H/h;
    const scalar eps2 = eps*eps;
    const scalar eps3 = eps*eps2;

    // const scalar C =
    //     sqrt(mag(g_)*h)*sqrt(1.0 + eps - 0.05*eps2 - (3.0/70.0)*eps3);
    const scalar C = this->celerity(H, h);

    // const scalar ts = 3.5*h/sqrt(eps);
    const scalar ts = this->waveLength(H, h) / 2.0;
    const scalar xa = -C*t + ts - X0 + x*cos(theta) + y*sin(theta);
    const scalar alfa = this->alfa(H, h);

    const scalar S = (1.0)/(cosh(alfa*(xa/h)));
    const scalar S2 = S*S;
    const scalar S4 = S2*S2;
    const scalar S6 = S2*S4;
    const scalar T = tanh(alfa*(xa/h));

    const scalar zbyh = z/h;
    const scalar zbyh2 = zbyh*zbyh;
    const scalar zbyh4 = zbyh2*zbyh2;
    /*
    scalar outa = eps*s2 - eps2*(-0.25*s2 + s4 + zbyh2*(1.5*s2 - 2.25*s4));
    scalar outb = 0.475*s2 + 0.2*s4 - 1.2*s6;
    scalar outc = zbyh2*(-1.5*s2 - 3.75*s4 + 7.5*s6);
    scalar outd = zbyh4*(-0.375*s2 + (45.0/16.0)*s4 - (45.0/16.0)*s6);
    
    scalar u = sqrt(mag(g_)*h)*(outa - eps3*(outb + outc + outd));
    */
    scalar ua = S2*eps - eps2*(-1.0/4.0*S2 + S4 + zbyh2*(3.0/2.0*S2-9.0/4.0*S4));
    scalar ub = -(19.0/40.0*S2 + 1.0/5.0*S4 - 6.0/5.0*S6);
    scalar uc = - zbyh2*( -3.0/2.0*S2-15.0/4.0*S4+15.0/2.0*S6 );
    scalar ud = - zbyh4*( -3.0/8.0*S2+45.0/16.0*S4-45.0/16.0*S6 ) ;

    scalar u = sqrt(mag(g_)*h) * (ua + eps3 * (ub + uc + ud));

    // decompose u into u and v
    const scalar v = u*sin(waveAngle_);
    u *= cos(waveAngle_);

    /*
    outa = eps*s2 - eps2*(0.375*s2 + 2*s4 + zbyh2*(0.5*s2 - 1.5*s4));
    outb = (49.0/640.0)*s2 - 0.85*s4 - 3.6*s6;
    outc = zbyh2*((-13.0/16.0)*s2 -(25.0/16.0)*s4 + 7.5*s6);
    outd = zbyh4*(-0.075*s2 -1.125*s4 - (27.0/16.0)*s6);
    // bug to calculate velocity w.
    const scalar w = sqrt(mag(g_)*h)*(outa - eps3*(outb + outc + outd));
    */
    scalar wa = - eps*S2 + eps2 * ( 3.0/8.0*S2 + 2.0*S4 + zbyh2*(1.0/2.0*S2 - 3.0/2.0*S4 ) );
    scalar wb = 49.0/640.0*S2 - 17.0/20.0*S4 - 18.0/5.0*S6 ;
    scalar wc = + zbyh2*(-13.0/16.0*S2 - 25.0/16.0*S4 + 15.0/2.0*S6 );
    scalar wd = + zbyh4*(-3.0/40.0*S2+9.0/8.0*S4-27.0/16.0*S6 ) ;

    scalar w = -sqrt(mag(g_)*h) * (sqrt(3.0*eps)*(zbyh)*T) * (wa + eps3 * (wb + wc + wd));

    return vector(u, v, w);
}


void Foam::waveModels::GrimshawNew::setLevel
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

Foam::waveModels::GrimshawNew::GrimshawNew
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
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::GrimshawNew::readDict(const dictionary& overrideDict)
{
    if (solitaryWaveModel::readDict(overrideDict))
    {
        return true;
    }

    return false;
}


void Foam::waveModels::GrimshawNew::setVelocity
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


void Foam::waveModels::GrimshawNew::info(Ostream& os) const
{
    solitaryWaveModel::info(os);
}


// ************************************************************************* //

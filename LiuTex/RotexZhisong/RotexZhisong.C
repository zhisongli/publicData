/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    stressComponents

Description
    Calculates and writes the scalar fields of the six components of the stress
    tensor sigma for each time.
    Zhisong @ SMU by 2023-12-08. lizhisongsjtu@163.com
    (0) This work is done based on the original Rotex.c and the work of Dongyue.
    (1) P, Q, R, delta are calculated with field manipulation functions. this saves the running time dramatically.
    (2) A scalar field, magRotex, is added to export the magnitude of Rotex.
    (3) Run: Rotex -latestTime | tee log.Rotex

    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
tensor rotation(vector u, vector v)
{
    scalar eps = 1.0e-10;
    tensor q(tensor::zero); 
    vector a(
               u.y()*v.z()-u.z()*v.y(),
               u.z()*v.x()-u.x()*v.z(),
               u.x()*v.y()-u.y()*v.x()
             );

   if(Foam::mag(a)<eps)
   {
         q.xx()= 1.0;
         q.yy()= 1.0;
         q.zz()= 1.0;
   }
   else
   {
         a=a/Foam::mag(a);
         scalar t=Foam::sign( u & v);
         scalar alpha = Foam::acos(t);
         scalar c = Foam::cos(alpha);
         scalar s = Foam::sin(alpha);

         q.xx()=Foam::pow(a.x(),2)*(1-c)+c;
         q.xy()=a.x()*a.y()*(1-c)-a.z()*s;
         q.xz()=a.x()*a.z()*(1-c)+a.y()*s;
 
         q.yx() = a.y()*a.x()*(1-c)+a.z()*s;
         q.yy() = Foam::pow(a.y(),2)*(1-c)+c;
         q.yz() = a.y()*a.z()*(1-c)-a.x()*s;

         q.zx() = a.z()*a.x()*(1-c)-a.y()*s;
         q.zy() = a.z()*a.y()*(1-c)+a.x()*s;
         q.zz() = Foam::pow(a.z(),2)*(1-c)+c;
   }
   return q;
}


int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );
        
        vector z0(0,0,1);

        // Check U exists
        if (Uheader.typeHeaderOk<volVectorField>(true))
        {
            mesh.readUpdate();

            Info<< "Reading U" << endl;
            volVectorField U(Uheader, mesh);
            volVectorField Rotex
            (
                IOobject
                (
                    "Rotex",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimless, vector::zero)
            );
            //
            volScalarField magRotex
            (
                IOobject
                (
                    "magRotex",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, scalar(0.0))
            );

            volTensorField gradU = fvc::grad(U);
            
            int myIndex = 0;
            forAll(U,cellI)
            {
                auto start = std::chrono::steady_clock::now();                                                   
                myIndex += 1;
                if (myIndex % 1000 == 0)
                {
                    Info << "myIndex = " << myIndex << endl;
                }
                // gradU
                tensor a=gradU[cellI];
                // tt is useless not
                // tensor tt =a & a;

                // tr(gradU), Eq4, P. use the default field manipulation function
                scalar aa1 = -tr(a);
                // scalar aa0=-(gradU.component(tensor::XX)()[cellI]
                //            + gradU.component(tensor::YY)()[cellI]
                //            + gradU.component(tensor::ZZ)()[cellI]
                //           );
                // Info << "aa0 - aa1 = " << aa0 - aa1 << endl;
                scalar aa = aa1;

                // Eq.5, Q, new expression from the function object Q.
                scalar bb1 = 0.5*(sqr(tr(a)) - tr(((a) & (a))));
                // scalar bb0 = -0.5*(tt.xx()+tt.yy()+tt.zz()-pow((tt.xx()+tt.yy()+tt.zz()),2));
                // Info << "bb0 - bb1 = " << bb0 - bb1 << endl;
                scalar bb = bb1;

                // Eq.6, R, 
                scalar cc1 = -det(a);
                // scalar cc0 = -( a.xx()*(a.yy()*a.zz()-a.yz()*a.zy())
                //               -a.xy()*(a.yx()*a.zz()-a.yz()*a.zx())
                //               +a.xz()*(a.yx()*a.zy()-a.yy()*a.zx())
                //              );
                // Info << "cc - cc1 = " << cc - cc1 << endl;
                scalar cc = cc1;

                // Eq.12, Qtidal
                // this typo is effectless for incompressible flow.
                // scalar qq=(Foam::pow(aa,2),3*bb)/9.0; 
                // Modifed to agree with the original expression
                scalar qq=(3*bb - Foam::pow(aa,2))/9.0;
                
                // Eq.12, Rtidal
                scalar rr=(2.0*Foam::pow(aa,3)-9.0*aa*bb+27.0*cc)/54.0;
                // Eq.12, calculate delta with the original expression, for simplicy
                // with the original expression, delta-delta1 is of order (-20), but Rx and Ry are changed
                scalar delta1 = Foam::pow(qq/3,3) + Foam::pow(rr/2,2);
                // scalar delta = 
                //     (
                //     18.0*aa*bb*cc - 
                //     4.0*Foam::pow(aa,3)*cc + 
                //     Foam::pow(aa,2)*pow(bb,2) -
                //     4.0*Foam::pow(bb,3) - 
                //     27.0*Foam::pow(cc,2)
                //     );
                // delta=-delta/108.0;
                scalar delta = delta1;
                // Info << "delta - delta1 = " << delta - delta1 << endl;

              
                if(delta > 0.0)
                {
                    scalar aaaa=-sign(rr)*Foam::pow(abs(rr)+Foam::sqrt(delta),scalar(1.0/3.0));
                    scalar bbbb=0;
 
                    if( aaaa !=0.0)
                       bbbb=qq/aaaa;
 
		            // scalar eig1c_r=-0.5*(aaaa+bbbb)-aa/3.0;
		            // scalar eig1c_i=0.5*Foam::sqrt(scalar(3.0))*(aaaa-bbbb);
		            // scalar eig2c_r=-0.5*(aaaa+bbbb)-aa/3.0;
		            // scalar eig2c_i=-0.5*Foam::sqrt(scalar(3.0))*(aaaa-bbbb);
                    scalar eig3r = aaaa+bbbb-aa/3.0;
                    scalar delta1(0.0),delta2(0.0),delta3(0.0);
                    delta1 = (a.xx()-eig3r)*(a.yy()-eig3r) -a.yx()*a.xy();
                    delta2 = (a.yy()-eig3r)*(a.zz()-eig3r) -a.yz()*a.zy();
                    delta3 = (a.xx()-eig3r)*(a.zz()-eig3r) -a.zx()*a.xz();
                  
                    if(delta1==0.0 && delta2==0.0 && delta3==0.0)
                    {
                        FatalErrorInFunction
	                      << "delta1-3 are:" << delta1 << ","
        	              << delta2 <<","<< delta3
                         << exit(FatalError);
                    }
                    vector vr(0,0,0);
                    
                    if(    Foam::mag(delta1)>=Foam::mag(delta2) 
                        &&  Foam::mag(delta1)>=Foam::mag(delta3)
                       ) 
                    {
                          vr.x()=(-(a.yy()-eig3r)*a.xz()+a.xy()*a.yz())/delta1;
                          vr.y()= (a.yx()*a.xz() - (a.xx()-eig3r)*a.yz())/delta1;
                          vr.z()=1.0;
                    }
                    else if(    Foam::mag(delta2)>=Foam::mag(delta1) 
                           &&  Foam::mag(delta2)>=Foam::mag(delta3)
                         )
                    {
                          vr.x()=1.0;
                          vr.y()=(-(a.zz()-eig3r)*a.yx()+a.yz()*a.zx())/delta2;
                          vr.z()=(a.zy()*a.yx() - (a.yy()-eig3r)*a.zx())/delta2;
                    }   
                    else if(    Foam::mag(delta3)>=Foam::mag(delta1) 
                           &&  Foam::mag(delta3)>=Foam::mag(delta2)
                         )
                    {
                          vr.x()=(-(a.zz()-eig3r)*a.xy()+a.xz()*a.zy())/delta3;
                          vr.y()= 1.0 ;
                          vr.z()=(a.zx()*a.xy() - (a.xx()-eig3r)*a.zy())/delta3;
                    }   
                    else
                         FatalErrorInFunction<< "vr error"<< exit(FatalError);
                 
                    vr = vr/Foam::sqrt(vr & vr);
 
                    tensor qqq=rotation(z0,vr);
                    tensor vg(qqq.T() & a);
                       vg=vg  & qqq; 
                              
                    scalar alpha = 0.5*Foam::sqrt(Foam::pow(vg.yy()-vg.xx(),2)+Foam::pow(vg.yx()+vg.xy(),2)); 
                    scalar beta = 0.5*(vg.yx()-vg.xy());
                 
                    if(Foam::magSqr(beta) > Foam::magSqr(alpha))
                    {  
                         scalar rm=0.0;
                         if(beta > 0.0)
                            rm=2.0*(beta-alpha);
                         else
                            rm=2.0*(beta+alpha);
                            
                         Rotex[cellI]=rm*vr; 
                         magRotex[cellI] = mag(Rotex[cellI]);

                    }
                }
                
                auto end = std::chrono::steady_clock::now();                                                     
                auto diff = end - start;

                // 
                //Info<< "Calculation using " 
                //    << std::chrono::duration <double, std::milli> (diff).count()                                 
                //    << " ms" << endl;
            }
            
            Rotex.write();
            magRotex.write();
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;

        }
        else
        {
            Info<< "    No U" << endl;
        }

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //

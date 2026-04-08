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
    Liutex

Description
    Calculates the Liutex vector field (third-generation vortex identification)
    based on the velocity gradient tensor.
    
        Zhisong @ SMU by 2023-12-08. lizhisongsjtu@163.com
        (0) This work is done based on the original Rotex.c and the work of Dongyue.
        (1) P, Q, R, delta are calculated with field manipulation functions. this saves the running time dramatically.
        (2) A scalar field, magRotex, is added to export the magnitude of Rotex.
        (3) Run: Rotex -latestTime | tee log.Rotex

    

    Liutex is a vector whose direction is the local rotation axis and whose
    magnitude is twice the local rigid-body angular speed.

    References:
      [1] Liu C. "Liutex - third generation of vortex definition and
          identification methods", Acta Aerodynamica Sinica, 2020, 38(3):
          413-431.  (Chinese review with full derivation)
      [2] Liu C., Yu Y. "Mathematical foundation of Liutex theory",
          J. Hydrodynamics, 2022, 34(6): 981-993.
      [3] Gao Y., Liu C. "Rortex and comparison with eigenvalue-based
          vortex identification criteria", Phys. Fluids, 2018, 30: 085107.
          (Explicit formula Eq. 35)
    

    Version history:
      RotexOrig.c      - original with 3 bugs (sign, Q-invariant, comma)
      RotexZhisong.C   - partial fix (Q/R via built-in, comma fix, sign unfixed)
      RotexZhisongV1.C - updated with all 3 bugs fixed, explicit Liutex formula, no Q-rotation et al.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

        if (Uheader.typeHeaderOk<volVectorField>(true))
        {
            mesh.readUpdate();
            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            // --- Output fields -------------------------------------------

            volVectorField Liutex
            (
                IOobject
                (
                    "Liutex",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedVector("zero", dimless, vector::zero)
            );

            volScalarField magLiutex
            (
                IOobject
                (
                    "LiutexMag",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, scalar(0.0))
            );

            // =============================================================
            //  Stage 1: Field-level computation of invariants
            //
            //  Velocity gradient tensor (OpenFOAM convention):
            //    gradU[i][j] = du_i / dx_j
            //  This is (∇v)^T in Liu's notation, which is the correct
            //  velocity gradient tensor matrix per Ref.[2] Definition 1.
            //
            //  Characteristic equation of gradU ([1] Eq. 3):
            //    λ³ + P·λ² + Q·λ + R = 0
            //
            //  The three invariants ([1] Eq. 4-6):
            //    P = -tr(A)
            //    Q = ½[tr(A)² - tr(A²)]
            //    R = -det(A)
            //
            //  Cardano's depressed cubic substitution λ = t - P/3:
            //    t³ + q'·t + r' = 0
            //  where ([1] near Eq. 12):
            //    q = (P² - 3Q) / 9
            //    r = (2P³ - 9PQ + 27R) / 54
            //    δ = r² - q³
            //
            //  δ > 0 ⟹ one real + two complex conjugate eigenvalues
            //        ⟹ local fluid rotation exists (vortex point)
            //
            //  All of the above are computed as field operations for
            //  vectorized memory access and compiler auto-vectorization,
            //  avoiding per-cell overhead of the forAll loop.
            // =============================================================

            volTensorField gradU(fvc::grad(U));

            // Cache A² = gradU · gradU (used in Q; avoids recomputation)
            volTensorField A2(gradU & gradU);

            volScalarField trA(tr(gradU));
            volScalarField fldP(-trA);                                // P = -tr(A)
            volScalarField fldQ(0.5*(sqr(trA) - tr(A2)));            // Q
            volScalarField fldR(-det(gradU));                         // R

            // Cardano intermediate quantities
            volScalarField fldqq((sqr(fldP) - 3.0*fldQ) / 9.0);     // q
            volScalarField fldrr                                      // r
            (
                (2.0*pow3(fldP) - 9.0*fldP*fldQ + 27.0*fldR) / 54.0
            );
            volScalarField fldDelta(sqr(fldrr) - pow3(fldqq));       // δ

            // Release intermediate tensor field to free memory
            A2.clear();

            // =============================================================
            //  Stage 2: Pre-filter vortex cells (δ > 0)
            //
            //  Only cells with δ > 0 have complex eigenvalues, meaning
            //  local fluid rotation exists. Non-vortex cells (δ ≤ 0)
            //  have three real eigenvalues and Liutex = 0 by definition.
            //
            //  In typical turbulent flows the vortex region is a small
            //  fraction of the total domain, so this filtering avoids
            //  expensive eigenvector computation on the majority of cells.
            // =============================================================

            // Estimate initial capacity: 10% of mesh as a heuristic
            DynamicList<label> vortexCells(mesh.nCells() / 10);

            const scalarField& deltaI = fldDelta.primitiveField();
            forAll(deltaI, cellI)
            {
                if (deltaI[cellI] > 0.0)
                {
                    vortexCells.append(cellI);
                }
            }

            Info<< "    Vortex cells: " << vortexCells.size()
                << " / " << mesh.nCells()
                << " (" << 100.0*vortexCells.size()/max(mesh.nCells(), 1)
                << "%)" << endl;

            // =============================================================
            //  Stage 3: Per-cell eigenvector and Liutex computation
            //
            //  Only executed on vortex cells identified in Stage 2.
            //  Pre-computed field values (qq, rr, delta) are reused
            //  to avoid redundant arithmetic.
            // =============================================================

            // Direct access to internal fields for efficiency
            const scalarField& qqI    = fldqq.primitiveField();
            const scalarField& rrI    = fldrr.primitiveField();
            const tensorField& gradUI = gradU.primitiveField();
            vectorField& LiutexI      = Liutex.primitiveFieldRef();
            scalarField& magLiutexI   = magLiutex.primitiveFieldRef();

            forAll(vortexCells, idx)
            {
                const label cellI = vortexCells[idx];
                const tensor& A   = gradUI[cellI];
                const scalar  qq  = qqI[cellI];
                const scalar  rr  = rrI[cellI];
                const scalar  del = deltaI[cellI]; // already > 0

                // ---------------------------------------------------------
                //  Cardano's formula for eigenvalues
                //
                //  AA = -sign(r)·(|r| + √δ)^{1/3}
                //  BB = q / AA
                //
                //  Real eigenvalue:
                //    λ_r = (AA + BB) - P/3
                //
                //  Complex conjugate pair:
                //    λ_cr = -½(AA + BB) - P/3
                //    λ_ci = ½√3·|AA - BB|
                // ---------------------------------------------------------

                const scalar sgnR = (rr >= 0.0) ? scalar(1) : scalar(-1);
                const scalar AA   = -sgnR * Foam::cbrt(mag(rr) + Foam::sqrt(del));
                const scalar BB   = (mag(AA) > SMALL) ? qq/AA : scalar(0);
                const scalar sumAB = AA + BB;

                // Real eigenvalue λ_r
                const scalar PP   = fldP[cellI];
                const scalar eigR = sumAB - PP/3.0;

                // Imaginary part of complex eigenvalue |λ_ci|
                const scalar lambda_ci = mag(0.5*Foam::sqrt(3.0)*(AA - BB));

                // ---------------------------------------------------------
                //  Real eigenvector: (A - λ_r·I)·r = 0
                //
                //  Solve via 2×2 minors. Three candidate 2×2 sub-determinants
                //  are formed by striking out each row/column in turn. The one
                //  with the largest absolute value is used to back-substitute
                //  for the eigenvector components, ensuring numerical stability.
                // ---------------------------------------------------------

                const scalar d1 = (A.xx() - eigR)*(A.yy() - eigR)
                                - A.yx()*A.xy();
                const scalar d2 = (A.yy() - eigR)*(A.zz() - eigR)
                                - A.yz()*A.zy();
                const scalar d3 = (A.xx() - eigR)*(A.zz() - eigR)
                                - A.zx()*A.xz();

                vector vr(vector::zero);
                const scalar ad1 = mag(d1);
                const scalar ad2 = mag(d2);
                const scalar ad3 = mag(d3);

                if (ad1 >= ad2 && ad1 >= ad3 && ad1 > SMALL)
                {
                    // Strike out row 3 (z-component free = 1)
                    vr.x() = (-(A.yy() - eigR)*A.xz() + A.xy()*A.yz()) / d1;
                    vr.y() = ( A.yx()*A.xz() - (A.xx() - eigR)*A.yz()) / d1;
                    vr.z() = 1.0;
                }
                else if (ad2 >= ad1 && ad2 >= ad3 && ad2 > SMALL)
                {
                    // Strike out row 1 (x-component free = 1)
                    vr.x() = 1.0;
                    vr.y() = (-(A.zz() - eigR)*A.yx() + A.yz()*A.zx()) / d2;
                    vr.z() = ( A.zy()*A.yx() - (A.yy() - eigR)*A.zx()) / d2;
                }
                else if (ad3 >= ad1 && ad3 >= ad2 && ad3 > SMALL)
                {
                    // Strike out row 2 (y-component free = 1)
                    vr.x() = (-(A.zz() - eigR)*A.xy() + A.xz()*A.zy()) / d3;
                    vr.y() = 1.0;
                    vr.z() = ( A.zx()*A.xy() - (A.xx() - eigR)*A.zy()) / d3;
                }
                else
                {
                    // All minors ≈ 0: degenerate (e.g., triple eigenvalue)
                    continue;
                }

                // Normalize to unit vector
                const scalar vrMag = mag(vr);
                if (vrMag < SMALL) continue;
                vr /= vrMag;

                // ---------------------------------------------------------
                //  Vorticity and direction convention
                //
                //  Vorticity ω = ∇ × v.  For A[i][j] = du_i/dx_j:
                //    ω_x = ∂w/∂y - ∂v/∂z = A.zy() - A.yz()
                //    ω_y = ∂u/∂z - ∂w/∂x = A.xz() - A.zx()
                //    ω_z = ∂v/∂x - ∂u/∂y = A.yx() - A.xy()
                //
                //  Convention ([1] below Eq. 38): ω·r > 0
                //  If not satisfied, flip r → -r.
                // ---------------------------------------------------------

                const vector omega
                (
                    A.zy() - A.yz(),
                    A.xz() - A.zx(),
                    A.yx() - A.xy()
                );

                scalar omegaDotR = omega & vr;
                if (omegaDotR < 0.0)
                {
                    vr = -vr;
                    omegaDotR = -omegaDotR;
                }

                // ---------------------------------------------------------
                //  Liutex magnitude ([1] Eq. 38, [3] Eq. 35):
                //
                //    R = (ω·r) - √[(ω·r)² - 4λ_ci²]
                //
                //  where r is the unit real eigenvector of (∇v)^T,
                //  ω is the vorticity, and λ_ci is the imaginary part
                //  of the complex eigenvalue pair.
                //
                //  The discriminant (ω·r)² - 4λ_ci² ≥ 0 is guaranteed
                //  analytically.  A negative value indicates numerical
                //  noise; in that case we fall back to R = ω·r.
                //
                //  Liutex vector:  R = R · r
                //  R/2 is the rigid-body angular speed of the local fluid.
                // ---------------------------------------------------------

                const scalar disc = sqr(omegaDotR) - 4.0*sqr(lambda_ci);
                const scalar Rmag = (disc >= 0.0)
                    ? omegaDotR - Foam::sqrt(disc)
                    : omegaDotR;

                LiutexI[cellI]    = Rmag * vr;
                magLiutexI[cellI] = mag(Rmag);
            }

            // Write output fields
            Liutex.write();
            magLiutex.write();

            Info<< "    ExecutionTime = " << runTime.elapsedCpuTime() << " s"
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

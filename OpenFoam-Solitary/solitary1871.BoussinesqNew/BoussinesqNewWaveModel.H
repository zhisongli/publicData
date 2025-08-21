/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 IH-Cantabria
    Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::waveModels::BoussinesqNew

Description
    BoussinesqNew wave model

\*---------------------------------------------------------------------------*/

#ifndef waveModels_BoussinesqNew_H
#define waveModels_BoussinesqNew_H

#include "solitaryWaveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{

/*---------------------------------------------------------------------------*\
                         Class BoussinesqNew Declaration
\*---------------------------------------------------------------------------*/

class BoussinesqNew
:
    public solitaryWaveModel
{
private:

    // Private Member Functions

        //- Wave height
        scalar eta
        (
            const scalar H,
            const scalar h,
            const scalar x,
            const scalar y,
            const scalar theta,
            const scalar t,
            const scalar X0
        ) const;

        //- Wave
        vector Deta
        (
            const scalar H,
            const scalar h,
            const scalar x,
            const scalar y,
            const scalar theta,
            const scalar t,
            const scalar X0
        ) const;

        //- Wave velocity
        vector Uf
        (
            const scalar H,
            const scalar h,
            const scalar x,
            const scalar y,
            const scalar theta,
            const scalar t,
            const scalar X0,
            const scalar z
        ) const;


protected:

    // Protected Member Functions

        //- Set the water level
        virtual void setLevel
        (
            const scalar t,
            const scalar tCoeff,
            scalarField& level
        ) const;

        //- Calculate the wave model velocity
        virtual void setVelocity
        (
            const scalar t,
            const scalar tCoeff,
            const scalarField& level
        );


public:

    //- Runtime type information
    TypeName("BoussinesqNew");

    //- Constructor
    BoussinesqNew
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const polyPatch& patch,
        const bool readFields = true
    );

    //- Destructor
    virtual ~BoussinesqNew() = default;


    // Public Member Functions

        //- Read from dictionary
        virtual bool readDict(const dictionary& overrideDict);

        //- Info
        virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

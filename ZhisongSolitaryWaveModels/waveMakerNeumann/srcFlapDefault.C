case motionTypes::flap:
{
    const pointField& points = patch().localPoints();
    scalarField motionX(patch().localPoints().size(), -1);

    forAll(points, pointi)
    {
        const label paddlei = pointToPaddle_[pointi];

        const scalar phaseTot =
            waveKx[paddlei]*xPaddle_[paddlei]
            + waveKy[paddlei]*yPaddle_[paddlei];

        const scalar depthRef = waterDepthRef_[paddlei];
        const scalar kh = waveK[paddlei]*depthRef;
        const scalar pz = points[pointi].component(2);

        const scalar m1 =
            (4*sinh(kh)/(sinh(2*kh) + 2*kh))
            * (sinh(kh) + 1/kh*(1 - cosh(kh)));

        const scalar boardStroke = waveHeight_/m1;

        motionX[pointi] = 0.5*boardStroke*sin(phaseTot - sigma*t);

        if (secondOrder_)
        {
            motionX[pointi] +=
                sqr(waveHeight_)/(16*depthRef)
                * (3*cosh(kh)/pow3(sinh(kh)) - 2/m1)
                * sin(phaseTot - 2*sigma*t);

        }

        motionX[pointi] *= 1.0 + (pz - zMinGb_ - depthRef)/depthRef;

    }

    Field<vector>::operator=(timeCoeff(t)*n_*motionX);

    break;
}

//
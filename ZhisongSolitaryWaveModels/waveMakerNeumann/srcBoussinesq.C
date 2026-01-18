case motionTypes::solitaryBoussinesq:
{
    const pointField& points = patch().localPoints();
    scalarField motionX(patch().localPoints().size(), -1);
    const scalar magG = mag(g());

    Info << "===== Boussinesq: t = " << t << nl;
    Info << "===== Boussinesq: order = " << order_ << nl;

    forAll(points, pointi)
    {
        const label paddlei = pointToPaddle_[pointi];
        const scalar depthRef = waterDepthRef_[paddlei];

        const scalar kappa = sqrt(0.75*waveHeight_/pow3(depthRef));
        const scalar celerity = sqrt(magG*(depthRef + waveHeight_));
        const scalar stroke = sqrt(16*waveHeight_*depthRef/3.0);
        const scalar hr = waveHeight_/depthRef;
        wavePeriod_ = 2.0/(kappa*celerity)*(3.8 + hr);
        const scalar tSolitary = -0.5*wavePeriod_ + t;

        // Newton-Raphson
        scalar theta1 = 0;
        scalar theta2 = 0;
        scalar er = 10000;
        const scalar error = 0.001;
        while (er > error)
        {
            theta2 =
                theta1
                - (theta1 - kappa*celerity*tSolitary + hr*tanh(theta1))
                /(1.0 + hr*(1.0/cosh(theta1))*(1.0/cosh(theta1)));

            er = mag(theta1 - theta2);
            theta1 = theta2;
        }

        motionX[pointi] =
            waveHeight_/(kappa*depthRef)*tanh(theta1) + 0.5*stroke;
    }

    Field<vector>::operator=(n_*motionX);

    break;
}

// 
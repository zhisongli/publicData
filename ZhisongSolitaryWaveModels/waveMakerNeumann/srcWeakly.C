case motionTypes::solitaryWeakly:
{
    const pointField& points = patch().localPoints();
    const scalar magG = mag(g());

    const scalar dtt = db().time().deltaT().value();

    const label order(order_); 
    const label StronglyNonlinear(0); // // 0-BoussinesqLaightoneGrimshawFenton, 1-RayleighChoi
    const label GR_MMT(Goring0MMT1_); // 0-GR, 1-MMT

    scalarField motionX(patch().localPoints().size(), -1);

    Info << "===================================================" << endl;

    Info << "waveHeight_ = " << waveHeight_ << nl;
    Info << "initialDepth_ = " << initialDepth_ << nl;

    Info << "# StronglyNonlinear = " << StronglyNonlinear << nl;
    Info << "    0: WeaklyNonlinear,   Boussinesq(1), Laiton(2), Grimshaw(3), Fenton(9), Qu(11) " << nl;
    Info << "    1: StronglyNonlinear, Rayleigh(1), Choi(2), Choi(3) " << nl;

    Info << "# order = " << order << nl;

    Info << nl;
    Info << "# Goring0MMT1 = " << Goring0MMT1_ << nl;
    Info << "    0: Goring method " << nl;
    Info << "    1: Malek-Mohammadi-Testik's " << nl;

    Info << "===================================================" << endl;
    // [1]Malek-Mohammadi S ,Testik Y F .New Methodology for Laboratory Generation of Solitary Waves[J].Journal of Waterway, Port, Coastal, and Ocean Engineering,2010,136(5):286-294. 

    
    forAll(points, pointi)
    {
        
        const label paddlei = pointToPaddle_[pointi];
        const scalar depthRef = waterDepthRef_[paddlei];

        const scalar epsilon(waveHeight_/depthRef);
        // calcualted with original equations
        const scalar kappa = _Kappa(waveHeight_, depthRef, order, StronglyNonlinear);
        const scalar celerity = _Celerity(waveHeight_, depthRef, order, StronglyNonlinear);
        // const scalar stroke = 2.0 * waveHeight_ / (kappa * depthRef);
        const scalar stroke = _Stroke(waveHeight_, depthRef, order, StronglyNonlinear);
        
        // wavePeriod_ = 2.0/(kappa*celerity)*(3.8 + epsilon);
        wavePeriod_ = _Period(waveHeight_, depthRef, order, StronglyNonlinear);

        const scalar tSolitary = -0.5*wavePeriod_ + t;

        const scalar tt = tSolitary - dtt;
        // const scalar tt = t - dtt;
        
        // Runge-Kutta
        scalar x = points[pointi].component(0) - xPaddle_[paddlei];

        scalar k1 = dtt * _dXidt(waveHeight_, depthRef, tt + 0.0 * dtt, x + 0.0*11, order, GR_MMT, StronglyNonlinear);
        scalar k2 = dtt * _dXidt(waveHeight_, depthRef, tt + 0.5 * dtt, x + 0.5*k1, order, GR_MMT, StronglyNonlinear);
        scalar k3 = dtt * _dXidt(waveHeight_, depthRef, tt + 0.5 * dtt, x + 0.5*k2, order, GR_MMT, StronglyNonlinear);
        scalar k4 = dtt * _dXidt(waveHeight_, depthRef, tt + 1.0 * dtt, x + 1.0*k3, order, GR_MMT, StronglyNonlinear);
        scalar xn1 = x + 1.0/6.0 *(k1 + 2.0*k2 + 2.0*k3 + k4);

        motionX[pointi] = xn1;

        
    }
    Info << "==========================================================" << nl;

    Field<vector>::operator=(n_*motionX);

    break;
}

//
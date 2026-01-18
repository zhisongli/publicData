
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_Kappa
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
)
{
    if (StronglyNonlinear == 0)
    {
        scalar k;
        scalar eps = H/h;
        scalar eps2(eps * eps);
        scalar eps3(eps2 * eps);

        const label N(10);
        scalar cKappa[N] = {1, -0.625, 0.554688, -0.561535, 0.567095, -0.602969, 0.624914, -0.670850, 0.700371, 0.0};
        scalar cK = 0;
        // for(label i=0; i < N; i++)
        for(label i=0; i < order; i++)
        {
            cK += cKappa[i] * pow(eps, i);
        }
        k = cK * sqrt(0.75*H/pow(h,3));

        // Info << "_Kappa: StronglyNonlinear = " << StronglyNonlinear << nl;

        return k;
    }
    else
    {
        vector tmpGammaEpsilon = _GammaEpsilon(H, h, order);
        scalar epsilon = tmpGammaEpsilon[1];
        scalar ks = sqrt(0.75*epsilon/(1.0+epsilon))/h;
        return ks;
    }
    
}

// 
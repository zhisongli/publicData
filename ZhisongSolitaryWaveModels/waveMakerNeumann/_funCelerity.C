
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_Celerity
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
)
{
    scalar C;
    const vector g_(g());
    if (StronglyNonlinear == 0)
        {
        scalar eps = H/h;
        scalar eps2(eps * eps);
        scalar eps3(eps2 * eps);
        
        const label N(10);
        scalar cC[N] = {1, 1, -0.05, -0.0428571, -0.0342857, -0.0315195, -0.0292784, -0.0268451, -0.0302634, -0.0219347};
        scalar cA = 0;
        for(label i=0; i < order+1; i++)
        {
            cA += cC[i] * pow(eps, i);
        }
        C = sqrt(cA) * sqrt(mag(g_)*h);
    }
    else
    {
        vector tmpGammaEpsilon = _GammaEpsilon(H, h, order);
        scalar gamma = tmpGammaEpsilon[0];
        scalar gamma1 = order >=2 ? gamma : 0;
        scalar gamma2 = order ==3 ? pow(gamma,2) : 0;
        scalar epsilon = tmpGammaEpsilon[1];
        scalar a = epsilon * h;
        scalar c0 = sqrt(mag(g_)*(a + h));
        // Info << "c0 = " << c0 << nl;
        C = 1 + 0.1 * epsilon * gamma1 + epsilon*(21*epsilon+40)*gamma2/1400;
        C *= c0;
    }
    
    return C;
}

//
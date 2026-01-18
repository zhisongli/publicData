
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_Eta
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  order,
    const label  StronglyNonlinear
)
{
    scalar eta(0);
    scalar eps = H/h;
    scalar eps2(eps * eps);
    scalar eps3(eps2 * eps);
    const vector g_(g());

    scalar C = _Celerity(H, h, order, StronglyNonlinear);
    scalar aux = _Kappa(H, h, order, StronglyNonlinear);
    scalar Xa = C*t - X;

    const scalar S = 1.0 / cosh(aux*Xa);
    const scalar T = tanh(aux*Xa);
    

    if (StronglyNonlinear == 0)
    {
        const label N(10), M(9);
        const scalar cEta[M][N] = 
        {
            {0, 1, -0.75, 0.625, -1.36817, 1.86057, -2.57413, 3.4572, -4.6849, 6.191},
            {0, 0, 0.75, -1.8875, 3.88033, -7.45136, 13.2856, -22.782, 37.670, -60.57},
            {0, 0, 0, 1.2625, -4.68304, 12.7637, -31.1191, 68.258, -139.28, 269.84},
            {0, 0, 0, 0, 2.17088, -11.4199, 40.1068, -116.974, 301.442, -712.125},
            {0, 0, 0, 0, 0, 4.24687, -28.4272, 120.49, -411.416, 1217.98},
            {0, 0, 0, 0, 0, 0, 8.728, -71.057, 355.069, -1384.37},
            {0, 0, 0, 0, 0, 0, 0, 18.608, -180.212, 1023.07},
            {0, 0, 0, 0, 0, 0, 0, 0, 41.412, -450.29},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 90.279}
        };

        for(label i=0; i < order+1; i++)
        {
            scalar tmpEta[N] = {0};
            for (label j = 0; j < i; j++)
            {
                tmpEta[i] += cEta[j][i] * pow(S, 2+2*j);
            }
            eta += tmpEta[i] * pow(eps, i);
        }
        eta *= h;
    }
    else
    {
        Info << "Strongly nonlinear model: eta ..." << nl;
        vector tmpGammaEpsilon = _GammaEpsilon(H, h, order);
        scalar gamma = tmpGammaEpsilon[0];
        scalar gamma1 = order >=2 ? gamma : 0;
        scalar gamma2 = order ==3 ? pow(gamma,2) : 0;
        scalar epsilon = tmpGammaEpsilon[1];
        scalar a = epsilon * h;
        scalar epsilonS(H/h);
        
        scalar a20 = 75 + 90*epsilon + 15*pow(epsilon,2);
        scalar a21 = -225 - 150*epsilon - 167*pow(epsilon,2);
        scalar a22 = -210*epsilon - 91*pow(epsilon,2);
        scalar a23 = 63*pow(epsilon,2);
        
        scalar a40 = 1653750 + 1527750*epsilon - 1794450*pow(epsilon,2) - 1935150*pow(epsilon,3) - 266700*pow(epsilon,4);
        scalar a41 = -18323550*epsilon - 9232110*pow(epsilon,2) - 11198450*pow(epsilon,3) - 1928556*pow(epsilon,4);
        scalar a42 = -12403125 - 21895650*epsilon - 24960330*pow(epsilon,2) + 15477950*pow(epsilon,3) + 3116512*pow(epsilon,4);
        scalar a43 = 16074450*epsilon - 75567510*pow(epsilon,2) - 78327450*pow(epsilon,3) - 24201716*pow(epsilon,4);
        scalar a44 = 146353500*pow(epsilon,2) + 108361950*pow(epsilon,3) + 26769680*pow(epsilon,4);
        scalar a45 = 9448950*pow(epsilon,3) + 8759170*pow(epsilon,4);
        scalar a46 = -4469535*pow(epsilon,4);

        scalar zeta0 = pow(S,2);
        scalar zeta2 = (gamma1*pow(S,2)*(a20 + (a21 + a22*pow(S,2) + a23*pow(S,4))*pow(T,2)))/300.;
        scalar zeta4 = (gamma2*pow(S,2)*(a40 + (a41 + pow(S,2)*(a42 + a43*pow(S,2) + a44*pow(S,4) + a45*pow(S,6) + a46*pow(S,8)))*pow(T,2)))/2.205e7;

        eta = (zeta0 + zeta2 + zeta4) * a;
    }

    
    return eta;
}

// 
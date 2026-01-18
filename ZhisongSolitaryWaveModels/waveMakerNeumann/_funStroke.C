
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_Stroke
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
)
{
    scalar stroke(0.0);
    const scalar kappa = _Kappa(H, h, order, StronglyNonlinear);
    const scalar celerity = _Celerity(H, h, order, StronglyNonlinear);
    
    if (StronglyNonlinear == 0)
    {
        const scalar eps(H / h);
        const label N(9);
        const scalar cStroke[N] = {1., -0.25, 0.04, -0.286502, 0.20549, -0.306553, 0.256789, -0.33643, 0.302132};

        for(label i=0; i < order; i++)
        {
            stroke += cStroke[i] * pow(eps, i);
        }
        stroke *= (2.0*eps/kappa);
    }
    else
    {
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
        
        scalar stroke1 = 1;
        scalar stroke2New = (a20 + a21/3. + (2*a22)/15. + (8*a23)/105.)/300 * gamma1;
        scalar stroke4New = (a40 + a41/3. + (2*a42)/15. + (8*a43)/105. + (16*a44)/315. + (128*a45)/3465. + (256*a46)/9009.)/22050000 * gamma2;

        stroke = (stroke1 + stroke2New + stroke4New) * (2.0*epsilon/kappa);
    }
   
    return stroke;
}

// 
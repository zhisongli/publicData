
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_Period
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
)
{
    scalar T(0.0);
    scalar eps = H/h;
    scalar eps2(eps * eps);
    scalar eps3(eps2 * eps);
    const vector g_(g());

    // wavePeriod_ = 2.0/(kappa*celerity)*(3.8 + hr);
    const scalar kappa = _Kappa(H, h, order, StronglyNonlinear);
    const scalar celerity = _Celerity(H, h, order, StronglyNonlinear);

    T = 2.0 / (kappa * celerity) * (3.8 + eps);

    return T;
}

// 
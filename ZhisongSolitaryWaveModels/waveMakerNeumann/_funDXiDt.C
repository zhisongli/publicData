
Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_dXidt
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  order,
    const label  GR_MMT,
    const label  StronglyNonlinear
)
{
    scalar eta(_Eta(H, h, t, X, order, StronglyNonlinear));
    scalar C(_Celerity(H, h, order, StronglyNonlinear));
    if (GR_MMT == 0)
    {
        return C * eta / (h + eta);
    }
    else if (GR_MMT == 1)
    {
        const vector g_(g());
        return sqrt(
            mag(g_) *
            eta / h *
            (h + 0.5 * eta) *
            eta / (h + eta)
        );
    }
}



Foam::scalar Foam::waveMakerNeumannPointPatchVectorField::_dXidtCF2003
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  GR_MMT,
    const Tensor<complex> k_nu_A_kk_C_nuk_u_uw
)
{
    const scalar x(X), y(0.0), theta(0.0), X0(0);

    const scalar tt(t / sqrt(h)); // convert from arbitary water depth h to h0=1.0
    scalar eta = etaCF2003(H/h, h/h, x/h, y/h, theta, tt, X0/h) * h;
    // convert celerity of h0=1.0 to the arbitary water depth.
    const scalar C = k_nu_A_kk_C_nuk_u_uw(1, 1).real() * sqrt(h);

    if (GR_MMT == 0)
    {
        return C * eta / (h + eta);
    }
    else if (GR_MMT == 1)
    {
        const vector g_(g());
        return sqrt(
            mag(g_) *
            eta / h *
            (h + 0.5 * eta) *
            eta / (h + eta)
        );
    }
}



// 

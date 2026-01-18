

// self defined variables
label order_;
label Goring0MMT1_;

// 
virtual vector _GammaEpsilon
(
    const scalar H,
    const scalar h,
    const label  order
);

//- 
virtual scalar _Kappa
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
);

virtual scalar _Celerity
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
);

virtual scalar _Period
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
);

virtual scalar _Stroke
(
    const scalar H,
    const scalar h,
    const label  order,
    const label  StronglyNonlinear
);

virtual scalar _Eta
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  order,
    const label  StronglyNonlinear
);


virtual scalar _dXidt
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  order,
    const label  GR_MMT,
    const label  StronglyNonlinear
);

virtual scalar _dXidtCF2003
(
    const scalar H,
    const scalar h,
    const scalar t,
    const scalar X,
    const label  GR_MMT,
    const Tensor<complex> k_nu_A_kk_C_nuk_u_uw
);


// 
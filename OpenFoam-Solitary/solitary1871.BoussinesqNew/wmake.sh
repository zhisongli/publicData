


rm -f ${FOAM_USER_LIBBIN}/libSolitaryBoussinesqABWaveModels.so


wmake libso

rm -rf lnInclude
rm -rf linux64GccDPInt32Opt

ls ${FOAM_USER_LIBBIN} -al
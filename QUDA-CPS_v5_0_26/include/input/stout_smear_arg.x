/*! Input paramters for AlgStoutSmear */
class StoutSmearArg
{     
  //! tolerance for the SU(3) projection
  Float tolerance;

  //! smear in hyper-plane orthoganal to this direction
  int   orthog;

  //! stout smearing coefficient 
  Float coef;
  
  memfun StoutSmearArg();
};

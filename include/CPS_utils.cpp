#include "CPS_utils.h"
#include <config.h>        //CPS Autoconf thing
#include <alg/wilson_matrix.h>

Float MMDag_re_tr(WilsonMatrix M) {

  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;

  Float re_tr = 0.0;

  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
        for(c2=0;c2<3;c2++){
          sc_idx = 2*(c2 + 3*s2 + 12*c1 + 36*s1);

          re_tr += M(s1,c1,s2,c2).real()*M(s1,c1,s2,c2).real() + M(s1,c1,s2,c2).imag()*M(s1,c1,s2,c2).imag();
        }
  return re_tr;
}

Float NMDag_re_tr(WilsonMatrix M, WilsonMatrix N) {

  int s1 = 0;
  int c1 = 0;
  int s2 = 0;
  int c2 = 0;
  int sc_idx = 0;

  Float re_tr = 0.0;

  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
        for(c2=0;c2<3;c2++){
          sc_idx = 2*(c2 + 3*s2 + 12*c1 + 36*s1);

          re_tr += M(s1,c1,s2,c2).real()*N(s2,c2,s1,c1).real() - M(s1,c1,s2,c2).imag()*N(s2,c2,s1,c1).imag();
        }
  return re_tr;
}


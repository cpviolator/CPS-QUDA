//CPS_utils.h
//generic utilities used by sum/ and 1pion/
#ifndef __CPS_UTILS__
#define __CPS_UTILS__

#include <config.h>        //CPS Autoconf thing
#include <alg/wilson_matrix.h>

USING_NAMESPACE_CPS

//Float MMDag_re_tr(WilsonMatrix M);
//Float NMDag_re_tr(WilsonMatrix M, WilsonMatrix N);

inline Float MMDag_re_tr(WilsonMatrix M) {
  int s1, c1, s2, c2;
  Float re_tr = 0.0;
  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
        for(c2=0;c2<3;c2++)
          re_tr += M(s1,c1,s2,c2).real()*M(s1,c1,s2,c2).real() + M(s1,c1,s2,c2).imag()*M(s1,c1,s2,c2).imag();
  return re_tr;
}

inline Float NMDag_re_tr(WilsonMatrix M, WilsonMatrix N) {
  int s1, c1, s2, c2;
  Float re_tr = 0.0;
  for(s1=0;s1<4;s1++)
    for(c1=0;c1<3;c1++)
      for(s2=0;s2<4;s2++)
        for(c2=0;c2<3;c2++)
          re_tr += M(s1,c1,s2,c2).real()*N(s2,c2,s1,c1).real() - M(s1,c1,s2,c2).imag()*N(s2,c2,s1,c1).imag();
  return re_tr;
}
#endif

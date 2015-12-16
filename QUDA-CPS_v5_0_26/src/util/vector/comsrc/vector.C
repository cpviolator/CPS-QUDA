#include<config.h>
CPS_START_NAMESPACE
/*!\file
  \brief  Definition of Vector and Matrix classes.

  Definitions of functions that perform operations on complex vectors.
  $Id: vector.C,v 1.12 2012-08-10 14:05:33 chulwoo Exp $
*/
//--------------------------------------------------------------------
//  CVS keywords
//
//  $Author: chulwoo $
//  $Date: 2012-08-10 14:05:33 $
//  $Header: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/vector.C,v 1.12 2012-08-10 14:05:33 chulwoo Exp $
//  $Id: vector.C,v 1.12 2012-08-10 14:05:33 chulwoo Exp $
//  $Name: not supported by cvs2svn $
//  $Locker:  $
//  $Revision: 1.12 $
//  $Source: /home/chulwoo/CPS/repo/CVS/cps_only/cps_pp/src/util/vector/comsrc/vector.C,v $
//  $State: Exp $
//
//--------------------------------------------------------------------

CPS_END_NAMESPACE
#include <util/vector.h>
#include <util/gjp.h>
#include <comms/glb.h>
CPS_START_NAMESPACE

//------------------------------------------------------------------
// 0, 1, 2
// 3, 4, 5
// 6, 7, 8
/*!
  \param m A linear array representation of a 3x3 complex matrix, such that 
  real part of the (i,j) element is at array position [6i+2j] 
  and the imaginary part of the (i,j) element is at array position [6i+2j+1].
  \post This matrix is the transpose of \a m.

  \a m must not be an alias of this matrix/
*/  
void Matrix::Trans(const IFloat *m) 
{
  Complex *dst = (Complex *)u;
  const Complex *src = (const Complex *)m;
  dst[0] = src[0]; dst[1] = src[3]; dst[2] = src[6];
  dst[3] = src[1]; dst[4] = src[4]; dst[5] = src[7];
  dst[6] = src[2]; dst[7] = src[5]; dst[8] = src[8];
}

/*!
  \return <em>|U^dagger U - I|^2</em>, where the norm used is the L2 norm
  of the matrix elements.
*/
IFloat Matrix::ErrorSU3() const
{
  Matrix tmp1, tmp2;

  tmp1.Dagger((const IFloat *)this);	// tmp1 = U~dag
  tmp2.DotMEqual(tmp1, *this); 
  
  tmp2 -= 1.0;	// tmp2 = U~dag U - I
  
  const IFloat *p = (const IFloat *)&tmp2;
  return dotProduct(p, p, COLORS*COLORS*2);
}

//Begin QUDA_CPS
// Use Cayley-Hamilton Theorem for SU(3) exp{iQ}.
// This algorithm is outlined in http://arxiv.org/pdf/hep-lat/0311018v1.pdf
// Equation numbers in the paper are referenced by [eq_no].
// src/util/vector/noarch/su3_util.C will come in handy too.

void Matrix::exp_iQ(void)
{  
  //Declarations
  Float inv3 = 1.0/3.0;
  Matrix Q;
  Q = *this;

  Float *p = (Float *)u;
  Float c0;
  Float c0_max;
  Float c1;
  Complex TRA;
  Float theta;
  Float u_p;  //u parameter
  Float w_p;  //w parameter
  Float *f0 = (Float*)smalloc(2*sizeof(Float));
  Float *f1 = (Float*)smalloc(2*sizeof(Float));
  Float *f2 = (Float*)smalloc(2*sizeof(Float));
  Matrix temp1;
  Matrix temp2;
  Float det[2];  

  //[14] c0 = det(Q) = 1/3Tr(Q^3)
  Q.Det(det);
  c0 = det[0];
  
  //[15] c1 = 1/2Tr(Q^2)
  // Q = Q^dag => Tr(Q^2) = Tr(QQ^dag) = sum_ab [Q_ab * Q_ab^*]
  temp1.DotMEqual(Q,Q);
  TRA = temp1.Tr();
  TRA *= 0.5;
  c1 = TRA.real();
  
  //We now have the coeffiecients c1 and c2.
  //We now find: exp(iQ) = f0*I + f1*Q + f2*Q^2
  //      where       fj = fj(c0,c1), j=0,1,2.
  
  //[17]
  c0_max  = 2*pow(c1*inv3,1.5);

  //[25]
  theta   = acos(c0/c0_max);

  //[23]
  u_p = sqrt(c1*inv3)*cos(theta*inv3);
  
  //[24]
  w_p = sqrt(c1)*sin(theta*inv3);
  
  //[29] Construct objects for fj = hj/(9u^2 - w^2).
  Float u_sq = u_p*u_p;
  Float w_sq = w_p*w_p;
  Float denom = 9*u_sq - w_sq;
  Float exp_iu_re = cos(u_p);
  Float exp_iu_im = sin(u_p);
  Float exp_2iu_re = exp_iu_re*exp_iu_re - exp_iu_im*exp_iu_im;
  Float exp_2iu_im = 2*exp_iu_re*exp_iu_im;
  Float cos_w = cos(w_p);
  Float sinc_w; sinc_w = Sinc(w_p);
  Float hj_re = 0.0;
  Float hj_im = 0.0;
  
  //[34] Test for c0 < 0.
  int parity = 0;
  if(c0 < 0) {
    c0 *= -1.0;
    parity = 1;
    //calculate fj with c0 > 0 and then convert all fj.
  }
  
  //Get all the numerators for fj,
  //[30] f0
  hj_re = (u_sq - w_sq)*exp_2iu_re + 8*u_sq*cos_w*exp_iu_re + 2*u_p*(3*u_sq + w_sq)*sinc_w*exp_iu_im;
  hj_im = (u_sq - w_sq)*exp_2iu_im - 8*u_sq*cos_w*exp_iu_im + 2*u_p*(3*u_sq + w_sq)*sinc_w*exp_iu_re;
  f0[0] = (hj_re/denom);
  f0[1] = (hj_im/denom);
  
  //[31] f1
  hj_re = 2*u_p*exp_2iu_re - 2*u_p*cos_w*exp_iu_re + (3*u_sq - w_sq)*sinc_w*exp_iu_im;
  hj_im = 2*u_p*exp_2iu_im + 2*u_p*cos_w*exp_iu_im + (3*u_sq - w_sq)*sinc_w*exp_iu_re;
  f1[0] = (hj_re/denom);
  f1[1] = (hj_im/denom);
  
  //[32] f2
  hj_re = exp_2iu_re - cos_w*exp_iu_re - 3*u_p*sinc_w*exp_iu_im;
  hj_im = exp_2iu_im + cos_w*exp_iu_im - 3*u_p*sinc_w*exp_iu_re;  
  f2[0] = (hj_re/denom);
  f2[1] = (hj_im/denom);
  
  //[34] If c0 < 0, apply tranformation  fj(-c0,c1) = (-1)^j f^*j(c0,c1)
  if (parity == 1) {
    f0[1] *= -1.0; 
    f1[0] *= -1.0;
    f2[1] *= -1.0;
  }
  
  //Create Complex data types for fj so we can use the Matrix
  //data type functions.
  Complex f0_c;
  Complex f1_c;
  Complex f2_c;
  
  f0_c.real(f0[0]);  
  f0_c.imag(f0[1]);

  f1_c.real(f1[0]);  
  f1_c.imag(f1[1]);  

  f2_c.real(f2[0]);  
  f2_c.imag(f2[1]);  
  
  //[19] Construct exp{iQ}
  Matrix exp_iQ;
  exp_iQ.ZeroMatrix();
  Matrix UNIT_mat;
  UNIT_mat = f0_c;
  // +f0*I
  exp_iQ += UNIT_mat;
  temp1 = Q;
  UNIT_mat = f1_c;
  temp2.DotMEqual(UNIT_mat,temp1);
  // +f1*Q
  exp_iQ += temp2;
  temp1.DotMEqual(Q,Q);
  UNIT_mat = f2_c;
  temp2.DotMEqual(UNIT_mat,temp1);
  // +f2Q^2
  exp_iQ += temp2;
  *this = exp_iQ;
  
  //Check that exp{iQ} is in SU(3)
  //Float err = exp_iQ.ErrorSU3();
  //printf("EXP: ERROR=%.16f\n",err);
  
  //free pointers
  sfree(f0);
  sfree(f1);
  sfree(f2);
}

//[33] Added one more term to the series given in the paper.
//EVAN W's CODE
IFloat Matrix::Sinc(IFloat x)
{
  if (x < 0.05 && x > -0.05) {
    //Let x = x^2
    x *= x;

    //1 - 1/6 x^2 (1 - 1/20 x^2 (1 - 1/42 x^2(1 - 1/72*x^2)))
    return 1.0-1.0/6.0*x*(1-0.05*x*(1-1.0/42.0*x*(1-1.0/72.0*x)));
  }
  else
    return sin(x)/x;
}
//End QUDA_CPS

//------------------------------------------------------------------
    /*!
      \param i The row index,
      \param j The column index,
      \return  The (i,j) matrix element
    */
Complex& Matrix::operator()(int i, int j)
{ return ((Complex*)u)[i*COLORS+j]; }


//------------------------------------------------------------------
    /*!
      \param i The row index,
      \param j The column index,
      \return  The (i,j) matrix element
    */
const Complex& Matrix::operator()(int i, int j) const
{ return ((Complex*)u)[i*COLORS+j]; }


//------------------------------------------------------------------
Complex Matrix::Char6() const
{
  Complex tmp(reChar6((IFloat *)u), imChar6((IFloat *)u)) ;
  return tmp ;
}

//------------------------------------------------------------------
Complex Matrix::Char8() const
{
  Complex tmp(reChar8((IFloat *)u)) ;
  return tmp ;
}

//------------------------------------------------------------------
Complex Matrix::Char10() const
{
  Complex tmp(reChar10((IFloat *)u), imChar10((IFloat *)u)) ;
  return tmp ;
}


//------------------------------------------------------------------
// The Vector class.
//------------------------------------------------------------------

//------------------------------------------------------------------
/*!
  \param len The number of real numbers in the vectors.
  \return The square norm of this vector summed over all nodes.
*/
//------------------------------------------------------------------
Float Vector::NormSqGlbSum(int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)&v, len);
  glb_sum_five((Float *)&sum);
  return Float(sum);
}
Float Vector::NormSqGlbSum4D(int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)&v, len);
  glb_sum((Float *)&sum);
  return Float(sum);
}

//------------------------------------------------------------------
/*!
  \param len The number of real numbers in the vectors. 
  For debugging purposes
*/
//-------------------------------------------------------------------

void Vector::Print(int len){
  for ( int n = 0; n < len; n+=2 ){
    printf("%dth vector, %e %e\n",n/2, v[n],v[n+1]);
  }

}
//------------------------------------------------------------------
/*!
  \param b Another vector
  \param len The number of real numbers in the vectors.
  \return The real part of the dot product (v,b) summed over all nodes.
*/
//------------------------------------------------------------------
Float Vector::ReDotProductGlbSum(const Vector *b, int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)b, len);
  glb_sum_five((Float *)&sum);
  return Float(sum);
}
Float Vector::ReDotProductGlbSum4D(const Vector *b, int len)
{
  IFloat sum = dotProduct((IFloat *)&v, (IFloat *)b, len);
  glb_sum((Float *)&sum);
  return Float(sum);
}

//------------------------------------------------------------------
/*!
  \param b Another vector
  \param len The number of real numbers in the vectors.
  \return The dot product of this vector with b, summed over all nodes.
 */
//------------------------------------------------------------------
Complex Vector::CompDotProductGlbSum(const Vector *b, int len)
{
  IFloat c_r, c_i;
  compDotProduct(&c_r, &c_i, (IFloat *)&v, (IFloat *)b, len);
  glb_sum_five((Float *)&c_r);
  glb_sum_five((Float *)&c_i);
  return Complex(c_r,c_i);
}
Complex Vector::CompDotProductGlbSum4D(const Vector *b, int len)
{
  IFloat c_r, c_i;
  compDotProduct(&c_r, &c_i, (IFloat *)&v, (IFloat *)b, len);
  glb_sum((Float *)&c_r);
  glb_sum((Float *)&c_i);
  return Complex(c_r,c_i);
}

//------------------------------------------------------------------
/*!
  For this array of vectors, with one on each lattice site, and
  some direction \a dir, sum the square norms of the vectors over the entire
  lattice on sites which
  have the same coordinate in direction \a dir. In other words, divide the
  global lattice into 3-dim slices perpendicular to direction \a dir and
  sum the square norms of the vectors in each slice.
  \param f_out The array holding the global sum for each slice.
  \param size The size of this vector on each lattice site.	
  \param dir The direction perpendicular to the slices.
*/
void Vector::NormSqArraySliceSum(Float *f_out, const int size, const int dir)
{
  char *cname = "Vector";
  char *fname = "NormSqArraySliceSum";
  Float *v1 = (Float *)smalloc(GJP.VolNodeSites()*sizeof(Float));
  if (v1 == 0)
    ERR.Pointer(cname, fname, "v1");
  VRB.Smalloc(cname, fname, "v1", v1, GJP.VolNodeSites()*sizeof(Float));

  IFloat *vv = (IFloat *)&v;
  for(int j=0; j < GJP.VolNodeSites(); j++, vv += size)
    v1[j] = Float(dotProduct(vv, vv, size));

  SliceArraySum(f_out, v1, dir);

  VRB.Sfree(cname, fname, "v1", v1);
  sfree(v1);
}

//------------------------------------------------------------------
/*!
  For an array of floating point numbers, one on each lattice site, and
  some direction \a dir, sum all the numbers over the entire lattice which
  have the same coordinate in direction \a dir. In other words, divide the
  global lattice into 3-dim slices perpendicular to direction \a dir and
  sum the array in each slice.
  \param sum The array holding the global sum for each slice.
  \param f_in The array to be sliced and summed
  \param dir The direction perpendicular to the slices.
*/
void Vector::SliceArraySum(Float *sum, const Float *f_in, const int dir)
{
  char *cname = "Vector";
  char *fname = "SliceArraySum";

  // Slice-sum the vector density to make a 1D vector
  if (dir < 0 || dir >= 4)
    ERR.General(cname,fname,"direction out of bounds");

  int nxx[4] = {GJP.Xnodes()*GJP.XnodeSites(),GJP.Ynodes()*GJP.YnodeSites(),
		GJP.Znodes()*GJP.ZnodeSites(),GJP.Tnodes()*GJP.TnodeSites()};
  int nx[4] = {GJP.XnodeSites(),GJP.YnodeSites(),
	       GJP.ZnodeSites(),GJP.TnodeSites()};
  int coord[4] = {GJP.XnodeCoor(),GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),GJP.TnodeCoor()};
  int len2 = 1;
  for(int i=0; i < dir; ++i)
    len2 *= nx[i];

  int len1 = len2 * nx[dir];
  int plane_offset = nx[dir]*coord[dir];

  int j;
  for(j=0; j < nxx[dir]; ++j)
    sum[j] = 0.0;

  for(j=0; j < GJP.VolNodeSites(); ++j)
  {
    int hypsec = j % len1;
    int plane  = hypsec / len2;
    sum[plane+plane_offset] += f_in[j];
  }

  for(j=0; j < nxx[dir]; ++j)
  {
    glb_sum(&(sum[j]));
  }
}


//------------------------------------------------------------------
/*!
  For an array of floating point numbers, one on each 5-dim lattice site, and
  some direction \a dir, sum all the numbers over the entire lattice which
  have the same coordinate in direction \a dir. In other words, divide the
  global lattice into 4-dim slices perpendicular to direction \a dir and
  sum the array in each slice.
  \param sum The array holding the global sum for each slice.
  \param f_in The array to be sliced and summed
  \param dir The direction perpendicular to the slices.
*/
//--------------------------------------------------------------------------
void Vector::SliceArraySumFive(Float *sum, const Float *f_in, const int dir)
{
  char *cname = "Vector";
  char *fname = "SliceArraySumFive";

  // Slice-sum the vector density to make a 1D vector
  if (dir < 0 || dir >= 5)
    ERR.General(cname,fname,"direction out of bounds");

  int nxx[5] = {GJP.Xnodes()*GJP.XnodeSites(),GJP.Ynodes()*GJP.YnodeSites(),
		GJP.Znodes()*GJP.ZnodeSites(),GJP.Tnodes()*GJP.TnodeSites(),
                GJP.Snodes()*GJP.SnodeSites()};
  int nx[5] = {GJP.XnodeSites(),GJP.YnodeSites(),
	       GJP.ZnodeSites(),GJP.TnodeSites(),
               GJP.SnodeSites()};
  int coord[5] = {GJP.XnodeCoor(),GJP.YnodeCoor(),
		  GJP.ZnodeCoor(),GJP.TnodeCoor(),
                  GJP.SnodeCoor()};
  int len2 = 1;
  for(int i=0; i < dir; ++i)
    len2 *= nx[i];

  int len1 = len2 * nx[dir];
  int plane_offset = nx[dir]*coord[dir];

  int j;
  for(j=0; j < nxx[dir]; ++j)
    sum[j] = 0.0;

  for(j=0; j < GJP.VolNodeSites()*GJP.SnodeSites(); ++j)
  {
    int hypsec = j % len1;
    int plane  = hypsec / len2;
    sum[plane+plane_offset] += f_in[j];
  }

  for(j=0; j < nxx[dir]; ++j)
  {
    glb_sum(&(sum[j]));
  }
}


//
//  y  +=  U x
//
void uDotXPlus(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y    += *u      * *x     - *(u+1)  * *(x+1) + *(u+2)  * *(x+2)
             - *(u+3)  * *(x+3) + *(u+4)  * *(x+4) - *(u+5)  * *(x+5);
    *(y+1)+= *u      * *(x+1) + *(u+1)  * *x     + *(u+2)  * *(x+3)
             + *(u+3)  * *(x+2) + *(u+4)  * *(x+5) + *(u+5)  * *(x+4);
    *(y+2)+= *(u+6)  * *x     - *(u+7)  * *(x+1) + *(u+8)  * *(x+2)
             - *(u+9)  * *(x+3) + *(u+10) * *(x+4) - *(u+11) * *(x+5);
    *(y+3)+= *(u+6)  * *(x+1) + *(u+7)  * *x     + *(u+8)  * *(x+3)
             + *(u+9)  * *(x+2) + *(u+10) * *(x+5) + *(u+11) * *(x+4);
    *(y+4)+= *(u+12) * *x     - *(u+13) * *(x+1) + *(u+14) * *(x+2)
             - *(u+15) * *(x+3) + *(u+16) * *(x+4) - *(u+17) * *(x+5);
    *(y+5)+= *(u+12) * *(x+1) + *(u+13) * *x     + *(u+14) * *(x+3)
             + *(u+15) * *(x+2) + *(u+16) * *(x+5) + *(u+17) * *(x+4);
}

/*! The 3x3 complex matrix is assumed to be stored in a linear form
  where the real part of the (i,j) element is at vector position [6i+2j]
  and the imaginary part of the (i,j) element is at vector position [6i+2j+1].
  \param y The vector \a y
  \param u The matrix \a M
  \param x The complex 3-vector
  \post \a y is the vector <em>y-Mx</em>

  The vector \a y must not alias vector \a x or \a u
*/
void uDotXMinus(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y    -= *u      * *x     - *(u+1)  * *(x+1) + *(u+2)  * *(x+2)
	     - *(u+3)  * *(x+3) + *(u+4)  * *(x+4) - *(u+5)  * *(x+5);
    *(y+1)-= *u      * *(x+1) + *(u+1)  * *x     + *(u+2)  * *(x+3)
	     + *(u+3)  * *(x+2) + *(u+4)  * *(x+5) + *(u+5)  * *(x+4);
    *(y+2)-= *(u+6)  * *x     - *(u+7)  * *(x+1) + *(u+8)  * *(x+2)
	     - *(u+9)  * *(x+3) + *(u+10) * *(x+4) - *(u+11) * *(x+5);
    *(y+3)-= *(u+6)  * *(x+1) + *(u+7)  * *x     + *(u+8)  * *(x+3)
	     + *(u+9)  * *(x+2) + *(u+10) * *(x+5) + *(u+11) * *(x+4);
    *(y+4)-= *(u+12) * *x     - *(u+13) * *(x+1) + *(u+14) * *(x+2)
	     - *(u+15) * *(x+3) + *(u+16) * *(x+4) - *(u+17) * *(x+5);
    *(y+5)-= *(u+12) * *(x+1) + *(u+13) * *x     + *(u+14) * *(x+3)
	     + *(u+15) * *(x+2) + *(u+16) * *(x+5) + *(u+17) * *(x+4);
}

/*! The 3x3 complex matrix is assumed to be stored in a linear form
  where the real part of the (i,j) element is at vector position [6i+2j]
  and the imaginary part of the (i,j) element is at vector position [6i+2j+1].
  \param y The vector \a y
  \param u The matrix \a M
  \param x The complex 3-vector \a x
  \post \a y is the vector <em>y + M^dagger x</em>

  The vector \a y must not alias vector \a x or \a u
*/
void uDagDotXPlus(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y    += *u      * *x     + *(u+1)  * *(x+1) + *(u+6)  * *(x+2)
	     + *(u+7)  * *(x+3) + *(u+12) * *(x+4) + *(u+13) * *(x+5);
    *(y+1)+= *u      * *(x+1) - *(u+1)  * *x     + *(u+6)  * *(x+3)
	     - *(u+7)  * *(x+2) + *(u+12) * *(x+5) - *(u+13) * *(x+4);
    *(y+2)+= *(u+2)  * *x     + *(u+3)  * *(x+1) + *(u+8)  * *(x+2)
	     + *(u+9)  * *(x+3) + *(u+14) * *(x+4) + *(u+15) * *(x+5);
    *(y+3)+= *(u+2)  * *(x+1) - *(u+3)  * *x     + *(u+8)  * *(x+3)
	     - *(u+9)  * *(x+2) + *(u+14) * *(x+5) - *(u+15) * *(x+4);
    *(y+4)+= *(u+4)  * *x     + *(u+5)  * *(x+1) + *(u+10) * *(x+2)
	     + *(u+11) * *(x+3) + *(u+16) * *(x+4) + *(u+17) * *(x+5);
    *(y+5)+= *(u+4)  * *(x+1) - *(u+5)  * *x     + *(u+10) * *(x+3)
	     - *(u+11) * *(x+2) + *(u+16) * *(x+5) - *(u+17) * *(x+4);
}

/*! The 3x3 complex matrix is assumed to be stored in a linear form
  where the real part of the (i,j) element is at vector position [6i+2j]
  and the imaginary part of the (i,j) element is at vector position [6i+2j+1].
  \param y The vector <em>M^dagger x</em>
  \param u The matrix \a M
  \param x The complex 3-vector \a x

  The vector \a y must not alias vector \a x or \a u
*/
void uDagDotXEqual(IFloat* y, const IFloat* u, const IFloat* x)
{
    *y     =  *u      * *x     + *(u+1)  * *(x+1) + *(u+6)  * *(x+2)
	    + *(u+7)  * *(x+3) + *(u+12) * *(x+4) + *(u+13) * *(x+5);
    *(y+1) =  *u      * *(x+1) - *(u+1)  * *x     + *(u+6)  * *(x+3)
	    - *(u+7)  * *(x+2) + *(u+12) * *(x+5) - *(u+13) * *(x+4);
    *(y+2) =  *(u+2)  * *x     + *(u+3)  * *(x+1) + *(u+8)  * *(x+2)
	    + *(u+9)  * *(x+3) + *(u+14) * *(x+4) + *(u+15) * *(x+5);
    *(y+3) =  *(u+2)  * *(x+1) - *(u+3)  * *x     + *(u+8)  * *(x+3)
	    - *(u+9)  * *(x+2) + *(u+14) * *(x+5) - *(u+15) * *(x+4);
    *(y+4) =  *(u+4)  * *x     + *(u+5)  * *(x+1) + *(u+10) * *(x+2)
	    + *(u+11) * *(x+3) + *(u+16) * *(x+4) + *(u+17) * *(x+5);
    *(y+5) =  *(u+4)  * *(x+1) - *(u+5)  * *x     + *(u+10) * *(x+3)
	    - *(u+11) * *(x+2) + *(u+16) * *(x+5) - *(u+17) * *(x+4);
}

CPS_END_NAMESPACE

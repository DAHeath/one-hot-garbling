#ifndef STANDARD_AES_SBOX_H__
#define STANDARD_AES_SBOX_H__


#include "share_matrix.h"


// See Logic Minimization Techniques with Applications to Cryptology
template <Mode mode>
ShareMatrix<mode> standard_aes_sbox(const ShareMatrix<mode>& x) {
  const auto x7 = x[0];
  const auto x6 = x[1];
  const auto x5 = x[2];
  const auto x4 = x[3];
  const auto x3 = x[4];
  const auto x2 = x[5];
  const auto x1 = x[6];
  const auto x0 = x[7];

  // top linear transformation
  const auto t0  = x1 ^ x2;
  const auto y1  = t0 ^ x7;
  const auto y14 = x3 ^ x5;
  const auto y8  = x0 ^ x5;
  const auto y4  = y1 ^ x3;
  const auto y5  = y1 ^ x6;
  const auto y13 = x0 ^ x6;
  const auto y12 = y13 ^ y14;
  const auto t1  = x4 ^ y12;
  const auto y15 = t1 ^ x5;
  const auto y10 = y15 ^ t0;
  const auto y20 = t1 ^ x1;
  const auto y9  = x0 ^ x3;
  const auto y11 = y20 ^ y9;
  const auto y17 = y10 ^ y11;
  const auto y16 = t0 ^ y11;
  const auto y21 = y13 ^ y16;
  const auto y3  =  y5 ^ y8;
  const auto y19 = y10 ^ y8;
  const auto y18 = x0 ^ y16;
  const auto y2  = y1 ^ x0;
  const auto y6  = y15 ^ x7;
  const auto y7  = x7 ^ y11;


  // middle non-linear section
  const auto t2  = y12 & y15;
  const auto t5  = y4 & x7;
  const auto t8  = y5 & y1;
  const auto t10 = y2 & y7;
  const auto t7  = y13 & y16;
  const auto t11 = t10 ^ t7;
  const auto t13 = y14 & y17;
  const auto t12 = y9 & y11;
  const auto t14 = t13 ^ t12;
  const auto t3  = y3 & y6;
  const auto t4  = t3 ^ t2;
  const auto t17 = t4 ^ t14;
  const auto t15 = y8 & y10;
  const auto t16 = t15 ^ t12;
  const auto t20 = t11 ^ t16;
  const auto t9  = t8 ^ t7;
  const auto t19 = t9 ^ t14;
  const auto t23 = t19 ^ y21;
  const auto t6  = t5 ^ t2;
  const auto t18 = t6 ^ t16;
  const auto t21 = t17 ^ y20;
  const auto t24 = t20 ^ y18;
  const auto t22 = t18 ^ y19;
  const auto t25 = t21 ^ t22;
  const auto t26 = t21 & t23;
  const auto t27 = t24 ^ t26;
  const auto t28 = t25 & t27;
  const auto t31 = t22 ^ t26;
  const auto t30 = t23 ^ t24;
  const auto t32 = t31 & t30;
  const auto t33 = t32 ^ t24;
  const auto t34 = t23 ^ t33;
  const auto t35 = t27 ^ t33;
  const auto t36 = t24 & t35;
  const auto t37 = t36 ^ t34;
  const auto t29 = t28 ^ t22;
  const auto t38 = t27 ^ t36;
  const auto t39 = t29 & t38;
  const auto t40 = t25 ^ t39;
  const auto t41 = t40 ^ t37;
  const auto t44 = t33 ^ t37;
  const auto z1  = t37 & y6;
  const auto z4  = t40 & y1;
  const auto t42 = t29 ^ t33;
  const auto t45 = t42 ^ t41;
  const auto z7  = t45 & y17;
  const auto z10 = t37 & y3;
  const auto z13 = t40 & y5;
  const auto z16 = t45 & y14;
  const auto z2  = t33 & x7;
  const auto z5  = t29 & y7;
  const auto z8  = t41 & y10;
  const auto z11 = t33 & y4;
  const auto z14 = t29 & y2;
  const auto z17 = t41 & y8;
  const auto t43 = t29 ^ t40;
  const auto z0  = t44 & y15;
  const auto z3  = t43 & y16;
  const auto z6  = t42 & y11;
  const auto z9  = t44 & y12;
  const auto z12 = t43 & y13;
  const auto z15 = t42 & y9;


  // bottom linear transformation
  const auto t46 = z15 ^ z16;
  const auto t49 = z9 ^ z10;
  const auto t52 = z7 ^ z8;
  const auto t55 = z16 ^ z17;
  const auto t58 = z4 ^ t46;
  const auto t54 = z6 ^ z7;
  const auto t59 = z3 ^ t54;
  const auto t64 = z4 ^ t59;
  const auto t63 = t49 ^ t58;
  const auto s0  = t59 ^ t63;
  const auto t50 = z2 ^ z12;
  const auto t53 = z0 ^ z3;
  const auto t57 = t50 ^ t53;
  const auto t61 = z14 ^ t57;
  const auto t62 = t52 ^ t58;
  const auto t65 = t61 ^ t62;
  const auto t67 = t64 ^ t65;
  const auto t47 = z10 ^ z11;
  const auto s5  = t47 ^ t65;
  const auto t48 = z5 ^ z13;
  const auto t56 = z12 ^ t48;
  const auto s6  = ~(t56 ^ t62);
  const auto t66 = z1 ^ t63;
  const auto s3  = t53 ^ t66;
  const auto s1  = ~(t64 ^ s3);
  const auto t51 = z2 ^ z5;
  const auto t60 = t46 ^ t57;
  const auto s7  = ~(t48 ^ t60);
  const auto s4  = t51 ^ t66;
  const auto s2  = ~(t55 ^ t67);


  ShareMatrix<mode> out(8, 1);
  out[7] = s0;
  out[6] = s1;
  out[5] = s2;
  out[4] = s3;
  out[3] = s4;
  out[2] = s5;
  out[1] = s6;
  out[0] = s7;

  return out;
}

#endif

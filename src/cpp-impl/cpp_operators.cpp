#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_param_mod.h"

using CppParamMod::NX_BLOCK;
using CppParamMod::NY_BLOCK;

void div(const int &k, 
    double (&div_out)[NY_BLOCK][NX_BLOCK],
    const double (&ux)[NY_BLOCK][NX_BLOCK],
    const double (&uy)[NY_BLOCK][NX_BLOCK]) {

  using CppGrid::hts;
  using CppGrid::htw;
  using CppGrid::kmt;
  using CppGrid::tarea_r;
  using CppConstantMod::C0;
  using CppConstantMod::P5;

  const int bid = 0;
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      div_out[j][i] = C0;
    }
  }
  for (int j = 1; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK-1; ++i) {
      if (k <= kmt[bid][j][i] - 1) {
        div_out[j][i] = P5 * (
            (ux[j  ][i+1] + ux[j-1][i+1]) * htw[bid][j  ][i+1]
          - (ux[j  ][i  ] + ux[j-1][i  ]) * htw[bid][j  ][i  ]
          + (uy[j  ][i+1] + uy[j  ][i  ]) * hts[bid][j  ][i  ]
          - (uy[j-1][i+1] + uy[j-1][i  ]) * hts[bid][j-1][i  ])
              * tarea_r[bid][j][i];
      }
    }
  }
  return ;
}

void grad(const int &k,
          double (&gradx)[NY_BLOCK][NX_BLOCK], 
          double (&grady)[NY_BLOCK][NX_BLOCK],
    const double (&f)[NY_BLOCK][NX_BLOCK]) {

	using CppConstantMod::C0;
	using CppConstantMod::P5;
	using CppGrid::kmu;
	using CppGrid::dyur;
	using CppGrid::dxur;

  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      gradx[j][i] = C0;
      grady[j][i] = C0;
    }
  }

  const int bid = 0;
  for (int j = 0; j < NY_BLOCK - 1; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      if (k <= (kmu[bid][j][i] - 1)) {
        gradx[j][i] = dxur[bid][j][i] * P5 
            * (f[j+1][i  ] - f[j][i-1] 
             - f[j+1][i-1] + f[j][i  ]);

        grady[j][i] = dyur[bid][j][i] * P5 
            * (f[j+1][i  ] - f[j][i-1]
             + f[j+1][i-1] - f[j][i  ]);
      }//end if
    }
  }//end for
  return ;
}

void zcurl(const int &k, 
    double (&curl)[NY_BLOCK][NX_BLOCK], 
    const double (&ux)[NY_BLOCK][NX_BLOCK], 
    const double (&uy)[NY_BLOCK][NX_BLOCK]) {
  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::kmt;
  using CppConstantMod::C0;
  using CppConstantMod::P5;
  const int bid = 0;
  for (int j = 0; j < NY_BLOCK; ++j) {
    for (int i = 0; i < NX_BLOCK; ++i) {
      curl[j][i] = C0;
    }
  }
  for (int j = 1; j < NY_BLOCK; ++j) {
    for (int i = 1; i < NX_BLOCK; ++i) {
      if (k <= kmt[bid][j][i] - 1) {
        curl[j][i] = P5 * (
            uy[j  ][i  ] * dyu[bid][j  ][i  ]
          + uy[j-1][i  ] * dyu[bid][j-1][i  ]
          - uy[j  ][i-1] * dyu[bid][j  ][i-1]
          - uy[j-1][i-1] * dyu[bid][j-1][i-1]
          - ux[j  ][i  ] * dxu[bid][j  ][i  ]
          - ux[j  ][i-1] * dxu[bid][j  ][i-1]
          + ux[j-1][i  ] * dxu[bid][j-1][i  ]
          + ux[j-1][i-1] * dxu[bid][j-1][i-1]);
      }
    }
  }
  return ;
}

#endif // LICOM_ENABLE_FORTRAN
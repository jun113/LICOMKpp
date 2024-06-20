#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_grid.h"
#include "../head/cpp_tracer_mod.h"

#include <string>
#include <cstdio>

using CppParamMod::KM;
using CppParamMod::JMT;
using CppParamMod::IMT;

template<typename T> T myMax(const T &val1, const T &val2) {
  return val1 > val2 ? val1 : val2;
}

template<typename T, typename... Args> T myMax(const T &val, const Args &... arg) {
  T result = myMax(arg...);
  return myMax(val, result);
}

template<typename T> T myMin(const T &val1, const T &val2) {
    return val1 < val2 ? val1 : val2;
}

template<typename T, typename... Args> T myMin(const T &val, const Args &... arg) {
  T result = myMin(arg...);
  return myMin(val, result);
}

void advection_momentum(
    const double (&uuu)[KM][JMT][IMT],
    const double (&vvv)[KM][JMT][IMT],
    const double (&www)[KM][JMT][IMT],
    double (&adv_uu)[KM][JMT][IMT],
    double (&adv_vv)[KM][JMT][IMT],
    const int &iblock) {

  using CppConstantMod::C0;
  using CppConstantMod::P5;
  using CppConstantMod::P25;
  using CppPconstMod::odzp;
  using CppPconstMod::adv_momentum;
  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::hue;
  using CppGrid::hun;
  using CppGrid::uarea_r;

  for (int k = 0; k < KM; ++k) {
    for (int j = 0; j < JMT; ++j) {
      for (int i = 0; i < IMT; ++i) {
        adv_uu[k][j][i] = C0;
        adv_vv[k][j][i] = C0;
      }
    }
  }
  std::string str_adv_momentum(adv_momentum);
  double u_wface[KM][JMT][IMT];
  double v_sface[KM][JMT][IMT];
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j][i-1] + uuu[k][j  ][i]) *
              P25 * hue[iblock][j  ][i-1];
          v_sface[k][j][i] = (vvv[k][j][i  ] + vvv[k][j+1][i]) *
              P25 * hun[iblock][j+1][i  ];
        }
      }
    }
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT-1; ++j) {
        for (int i = 1; i < IMT; ++i) {
          u_wface[k][j][i] = 
              (uuu[k][j  ][i-1] * dyu[iblock][j  ][i-1] +
               uuu[k][j  ][i  ] * dyu[iblock][j  ][i  ]) * P25;
          v_sface[k][j][i] = 
              (vvv[k][j  ][i  ] * dxu[iblock][j  ][i  ] +
               vvv[k][j+1][i  ] * dxu[iblock][j+1][i  ]) * P25;
        }
      }
    }
  } else {
    if (CppParamMod::mytid == 0) {
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  double adv_z1, adv_z2, adv_z3, adv_z4;
  if (str_adv_momentum.find("centered") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_uu[k][j][i] = (
              - u_wface[k][j  ][i  ] * (uuu[k][j  ][i  ] - uuu[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (uuu[k][j  ][i+1] - uuu[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (uuu[k][j+1][i  ] - uuu[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (uuu[k][j  ][i  ] - uuu[k][j-1][i  ])) 
                  * uarea_r[iblock][j][i];

          adv_vv[k][j][i] = (
              - u_wface[k][j  ][i  ] * (vvv[k][j  ][i  ] - vvv[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (vvv[k][j  ][i+1] - vvv[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (vvv[k][j+1][i  ] - vvv[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (vvv[k][j  ][i  ] - vvv[k][j-1][i  ])) 
                  * uarea_r[iblock][j][i];

          if (k == 0) {
            adv_z1 = 0.0;
            adv_z3 = 0.0;
          } else {
            adv_z1 = www[k][j][i] * (uuu[k-1][j][i] - uuu[k][j][i]);
            adv_z3 = www[k][j][i] * (vvv[k-1][j][i] - vvv[k][j][i]);
          }

          if (k == KM-1) {
            adv_z2 = 0.0;
            adv_z4 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (uuu[k][j][i] - uuu[k+1][j][i]);
            adv_z4 = www[k+1][j][i] * (vvv[k][j][i] - vvv[k+1][j][i]);
          }

          adv_uu[k][j][i] -= P5 * odzp[k] *
              (adv_z1 + adv_z2);
          adv_vv[k][j][i] -= P5 * odzp[k] *
              (adv_z3 + adv_z4);
        }
      }
    }
  } else if (str_adv_momentum.find("flux") != str_adv_momentum.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_uu[k][j][i] = (
              - u_wface[k][j  ][i  ] * (uuu[k][j  ][i  ] + uuu[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (uuu[k][j  ][i+1] + uuu[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (uuu[k][j  ][i  ] + uuu[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (uuu[k][j+1][i  ] + uuu[k][j  ][i  ])) 
                  * uarea_r[iblock][j][i];

          adv_vv[k][j][i] = (
              - u_wface[k][j  ][i  ] * (vvv[k][j  ][i  ] + vvv[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (vvv[k][j  ][i+1] + vvv[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (vvv[k][j  ][i  ] + vvv[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (vvv[k][j+1][i  ] + vvv[k][j  ][i  ])) 
                  * uarea_r[iblock][j][i];
   
          if (k == 0) {
            adv_z1 = 0.0;
            adv_z3 = 0.0;
          } else {
            adv_z1 = www[k][j][i] * (
                uuu[k-1][j][i] + uuu[k][j][i]) * P5;
            adv_z3 = www[k][j][i] * (
                vvv[k-1][j][i] + vvv[k][j][i]) * P5;
          }
 
          if (k == KM-1) {
            adv_z2 = 0.0;
            adv_z4 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (
                uuu[k][j][i] + uuu[k+1][j][i]) * P5; 
            adv_z4 = www[k+1][j][i] * (
                vvv[k][j][i] + vvv[k+1][j][i]) * P5; 
          }
          adv_uu[k][j][i] -= odzp[k] *
              (adv_z2 - adv_z1);
          adv_vv[k][j][i] -= odzp[k] *
              (adv_z4 - adv_z3);
        }
      }
    }
  } else {
    if (CppParamMod::mytid == 0) {
      printf("The false advection option for tracer\n");
    }
    exit(0);
  }
  return ;
}

void advection_tracer(const double (&uuu)[KM][JMT][IMT],
    const double (&vvv)[KM][JMT][IMT], const double (&www)[KM+1][JMT][IMT],
        const double (&ttt)[KM][JMT][IMT], double (&adv_tt)[KM][JMT][IMT],
            const int &iblock, const int &mtracer) {
  using CppGrid::dxu;
  using CppGrid::dyu;
  using CppGrid::hue;
  using CppGrid::hun;
  using CppGrid::htw;
  using CppGrid::hts;
  using CppGrid::tarea_r;
  using CppPconstMod::dts;
  using CppPconstMod::vit;
  using CppPconstMod::nss;
  using CppPconstMod::odzp;
  using CppPconstMod::odzt;
  using CppPconstMod::adv_tracer;
  using CppTracerMod::ax;
  using CppTracerMod::ay;
  using CppTracerMod::az;
  using CppConstantMod::P5;
  using CppConstantMod::P25;

  double u_wface[KM][JMT][IMT], v_sface[KM][JMT][IMT];

  const std::string str_adv_tracer(adv_tracer);
  if (str_adv_tracer.find("centered") != str_adv_tracer.npos ||
      str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          v_sface[k][j][i] = (vvv[k][j][i] + vvv[k][j][i+1]) 
              * hts[iblock][j][i] * P25;
        }
      }
    }
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 0; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j-1][i] + uuu[k][j][i]) 
              * htw[iblock][j][i] * P25;
        }
      }
    }
  }

  if (str_adv_tracer.find("flux") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 0; j < JMT; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          v_sface[k][j][i] = (vvv[k][j][i  ] * dxu[iblock][j][i  ]
                            + vvv[k][j][i+1] * dxu[iblock][j][i+1]) * P25;
        }
      }
    }
    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 0; i < IMT; ++i) {
          u_wface[k][j][i] = (uuu[k][j-1][i] * dyu[iblock][j-1][i]
                            + uuu[k][j  ][i] * dyu[iblock][j  ][i]) * P25; 
        }
      }
    }
  }


  if (str_adv_tracer.find("centered") != str_adv_tracer.npos) {
    // k = 0
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {
        adv_tt[0][j][i] = (
            - u_wface[0][j  ][i  ] * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
            - u_wface[0][j  ][i+1] * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ])
            - v_sface[0][j  ][i  ] * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
            - v_sface[0][j-1][i  ] * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]))
                * tarea_r[iblock][j][i];
        
        ax[iblock][mtracer][0][j][i] += (
            - u_wface[0][j  ][i  ] * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
            - u_wface[0][j  ][i+1] * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ]))
                * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        ay[iblock][mtracer][0][j][i] += (
            - v_sface[0][j  ][i  ] * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
            - v_sface[0][j-1][i  ] * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]))
                * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        const double adv_z2 = www[1][j][i] * (ttt[0][j][i] - ttt[1][j][i]);

        adv_tt[0][j][i] -= P5 * odzp[0] * adv_z2;

        az[iblock][mtracer][0][j][i] -= P5 * odzp[0] * adv_z2 
            / static_cast<double>(nss); 
      }
    }

    for (int k = 1; k < KM-1; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_tt[k][j][i] = (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ])
              - v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]))
                  * tarea_r[iblock][j][i];
      
          ax[iblock][mtracer][k][j][i] += (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
              - u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]))
                  * tarea_r[iblock][j][i] / static_cast<double>(nss);

          ay[iblock][mtracer][k][j][i] += (
              - v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]))
                  * tarea_r[iblock][j][i] / static_cast<double>(nss);
               
          const double adv_z1 = www[k][j][i] 
              * (ttt[k-1][j][i] - ttt[k  ][j][i]);
               
          const double adv_z2 = www[k+1][j][i] 
              * (ttt[k  ][j][i] - ttt[k+1][j][i]);

          adv_tt[k][j][i] -= P5 * odzp[k] * (adv_z1 + adv_z2);

          az[iblock][mtracer][k][j][i] -= P5 * odzp[k]
              * (adv_z1 + adv_z2) / static_cast<double>(nss);
        }
      }
    }
    // k = KM - 1
    for (int j = 2; j < JMT-2; ++j) {
      for (int i = 2; i < IMT-2; ++i) {
        adv_tt[KM-1][j][i] = (
            - u_wface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
            - u_wface[KM-1][j  ][i+1] 
                * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j-1][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]))
                    * tarea_r[iblock][j][i];
        
        ax[iblock][mtracer][KM-1][j][i] += (
            - u_wface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
            - u_wface[KM-1][j  ][i+1] 
                * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ]))
                    * tarea_r[iblock][j][i] / static_cast<double>(nss);
                 
        ay[iblock][mtracer][KM-1][j][i] += (
            - v_sface[KM-1][j  ][i  ] 
                * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
            - v_sface[KM-1][j-1][i  ] 
                * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]))
                    * tarea_r[iblock][j][i] / static_cast<double>(nss);

        const double adv_z1 = www[KM-1][j][i] 
            * (ttt[KM-2][j][i] - ttt[KM-1][j][i]);

        adv_tt[KM-1][j][i] -= P5 * odzp[KM-1] * adv_z1;

        az[iblock][mtracer][KM-1][j][i] -= P5 * odzp[KM-1] * adv_z1
            / static_cast<double>(nss); 
      }
    }
  } else if (str_adv_tracer.find("flux") != str_adv_tracer.npos) {
    for (int k = 0; k < KM; ++k) {
      for (int j = 2; j < JMT-2; ++j) {
        for (int i = 2; i < IMT-2; ++i) {
          adv_tt[k][j][i] = (
              - u_wface[k][j  ][i  ] * (ttt[k][j  ][i  ] + ttt[k][j  ][i-1])
              + u_wface[k][j  ][i+1] * (ttt[k][j  ][i+1] + ttt[k][j  ][i  ])
              - v_sface[k][j-1][i  ] * (ttt[k][j  ][i  ] + ttt[k][j-1][i  ])
              + v_sface[k][j  ][i  ] * (ttt[k][j+1][i  ] + ttt[k][j  ][i  ]))
                  * tarea_r[iblock][j][i];
          double adv_z1, adv_z2; 
          if (k == 0) {
            adv_z1 = 0.0;
          } else {
            adv_z1 = www[k  ][j][i] * (ttt[k-1][j][i] + ttt[k  ][j][i]) * P5;
          }
          if (k == KM-1) {
            adv_z2 = 0.0;
          } else {
            adv_z2 = www[k+1][j][i] * (ttt[k  ][j][i] + ttt[k+1][j][i]) * P5;
          }
          adv_tt[k][j][i] -= odzp[k] * (adv_z2 - adv_z1);
        }
      }
    }
  } else if (str_adv_tracer.find("tspas") != str_adv_tracer.npos) {
    double at00[KM][JMT][IMT], atmax[KM][JMT][IMT], atmin[KM][JMT][IMT];

    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          const double adv_x0 = (
              ttt[k][j][i+1] + ttt[k][j][i  ]) 
                  * u_wface[k][j][i+1] * tarea_r[iblock][j][i]
           - (ttt[k][j][i  ] + ttt[k][j][i-1]) 
                  * u_wface[k][j][i  ] * tarea_r[iblock][j][i];


          const double temp_1 = (ttt[k][j+1][i] + ttt[k][j  ][i]) 
              * v_sface[k][j  ][i];
          const double temp_2 = (ttt[k][j  ][i] + ttt[k][j-1][i]) 
              * v_sface[k][j-1][i];
          const double adv_y0 = (temp_1 - temp_2) * tarea_r[iblock][j][i];
          // TODO wjl 20211116
          /*
          const double adv_y0 = (
              (ttt[k][j+1][i] + ttt[k][j  ][i]) * v_sface[k][j  ][i]
            - (ttt[k][j  ][i] + ttt[j][j-1][i]) * v_sface[k][j-1][i])
                * tarea_r[iblock][j][i];
          */

          const double adv_xy1 = - dts * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);

          const double adv_xy2 =   dts * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);

          const double adv_xy3 = - dts * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(v_sface[k][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);

          const double adv_xy4 =   dts * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]) 
              * 2.0 * tarea_r[iblock][j][i] * pow(v_sface[k][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
             
          const double adv_c1 = - ttt[k][j][i] 
              * (u_wface[k][j  ][i+1] - u_wface[k][j  ][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                   
          const double adv_c2 = - ttt[k][j][i] 
              * (v_sface[k][j  ][i  ] - v_sface[k][j-1][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                    
          double adv_za, adv_zc;
          double adv_zb1, adv_zb2;
          if (k == 0) {
            adv_za = - 0.5 * odzp[0] * www[1][j][i] 
                * (ttt[1][j][i] + ttt[0][j][i]);
                 
            adv_zb1 = 0.0;

            adv_zb2 = 0.5 * odzp[0] * pow(www[1][j][i], 2) * odzt[1]
                * (ttt[0][j][i] - ttt[1][j][i]) * dts;

            adv_zc = odzp[0] * ttt[0][j][i] * www[1][j][i];
          } else if (k == KM-1) {

            adv_za = 0.5 * odzp[KM-1] * www[KM-1][j][i]
                * (ttt[KM-1][j][i] + ttt[KM-2][j][i]);

            adv_zb1 = - 0.5 * odzp[KM-1] * pow(www[KM-1][j][i], 2) * odzt[KM-1]
                * (ttt[KM-2][j][i] - ttt[KM-1][j][i]) * dts;

            adv_zb2 = 0.0;

            adv_zc = - odzp[KM-1] * ttt[KM-1][j][i] * www[KM-1][j][i];
          } else {
            adv_za = 0.5 * odzp[k] * www[k  ][j][i] 
                       * (ttt[k][j][i] + ttt[k-1][j][i])
                   - 0.5 * odzp[k] * www[k+1][j][i] 
                       * (ttt[k][j][i] + ttt[k+1][j][i]);

            adv_zb1 = - 0.5 * odzp[k] * pow(www[k  ][j][i], 2) * odzt[k  ]
                * (ttt[k-1][j][i] - ttt[k  ][j][i]) * dts;

            adv_zb2 =   0.5 * odzp[k] * pow(www[k+1][j][i], 2) * odzt[k+1]
                * (ttt[k  ][j][i] - ttt[k+1][j][i]) * dts;
                
            adv_zc = - odzp[k] * ttt[k][j][i] 
                * (www[k  ][j][i] - www[k+1][j][i]);

          }
          const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
          const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
          const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

          at00[k][j][i] = ttt[k][j][i] + (adv_xx + adv_yy + adv_zz) * dts;
        }
      }
    }

    const double wt1 = - 1.0e10;
    const double wt2 = + 1.0e10;

    for (int k = 0; k < KM; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          if (k == 0) {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt2);
          } else if (k == KM-1) {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt2);
          } else {
            atmax[k][j][i] = myMax(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt1, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt1, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt1, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt1, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt1, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt1, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt1);

            atmin[k][j][i] = myMin(
                ttt[k  ][j  ][i  ] * vit[iblock][k  ][j  ][i  ] 
                    + (1.0 - vit[iblock][k  ][j  ][i  ]) * wt2, 
                ttt[k  ][j-1][i  ] * vit[iblock][k  ][j-1][i  ] 
                    + (1.0 - vit[iblock][k  ][j-1][i  ]) * wt2, 
                ttt[k  ][j+1][i  ] * vit[iblock][k  ][j+1][i  ] 
                    + (1.0 - vit[iblock][k  ][j+1][i  ]) * wt2, 
                ttt[k  ][j  ][i-1] * vit[iblock][k  ][j  ][i-1] 
                    + (1.0 - vit[iblock][k  ][j  ][i-1]) * wt2, 
                ttt[k  ][j  ][i+1] * vit[iblock][k  ][j  ][i+1] 
                    + (1.0 - vit[iblock][k  ][j  ][i+1]) * wt2, 
                ttt[k-1][j  ][i  ] * vit[iblock][k-1][j  ][i  ] 
                    + (1.0 - vit[iblock][k-1][j  ][i  ]) * wt2, 
                ttt[k+1][j  ][i  ] * vit[iblock][k+1][j  ][i  ] 
                    + (1.0 - vit[iblock][k+1][j  ][i  ]) * wt2);
          }
        }
      }
    }
    // k = 0
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        double adv_xy1;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j  ][i+1] > atmax[0][j  ][i+1] || 
            at00[0][j  ][i+1] < atmin[0][j  ][i+1]) {

          adv_xy1 = - (ttt[0][j  ][i+1] - ttt[0][j  ][i  ])
              * fabs(u_wface[0][j  ][i+1]) * tarea_r[iblock][j][i];
        } else {
          adv_xy1 = - dts * (ttt[0][j  ][i+1] - ttt[0][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[0][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
        }

        double adv_xy2;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j  ][i-1] > atmax[0][j  ][i-1] || 
            at00[0][j  ][i-1] < atmin[0][j  ][i-1]) {

          adv_xy2 =   (ttt[0][j  ][i  ] - ttt[0][j  ][i-1])
              * fabs(u_wface[0][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy2 =   dts * (ttt[0][j  ][i  ] - ttt[0][j  ][i-1]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[0][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
        }

        double adv_xy3;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j+1][i  ] > atmax[0][j+1][i  ] || 
            at00[0][j+1][i  ] < atmin[0][j+1][i  ]) {

          adv_xy3 = - (ttt[0][j+1][i  ] - ttt[0][j  ][i  ])
              * fabs(v_sface[0][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy3 = - dts * (ttt[0][j+1][i  ] - ttt[0][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[0][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
        }

        double adv_xy4;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[0][j-1][i  ] > atmax[0][j-1][i  ] || 
            at00[0][j-1][i  ] < atmin[0][j-1][i  ]) {
 
          adv_xy4 =   (ttt[0][j  ][i  ] - ttt[0][j-1][i  ])
              * fabs(v_sface[0][j-1][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy4 =   dts * (ttt[0][j  ][i  ] - ttt[0][j-1][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[0][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
        }

        const double adv_zb1 = 0.0;
        double adv_zb2;
        if (at00[0][j  ][i  ] > atmax[0][j  ][i  ] || 
            at00[0][j  ][i  ] < atmin[0][j  ][i  ] || 
            at00[1][j  ][i  ] > atmax[1][j  ][i  ] || 
            at00[1][j  ][i  ] < atmin[1][j  ][i  ]) {
 
          adv_zb2 = 0.5 * fabs(www[1][j][i]) * odzp[0]
              * (ttt[0][j][i] - ttt[1][j][i]);
        } else {
          adv_zb2 = 0.5 * odzp[0] * pow(www[1][j][i], 2)
              * odzt[1] * (ttt[0][j][i] - ttt[1][j][i]) * dts;
        }

        const double adv_za = - 0.5 * odzp[0] * www[1][j][i]
            * (ttt[1][j][i] + ttt[0][j][i]);
        const double adv_zc = odzp[0] * ttt[0][j][i] * www[1][j][i];

        const double adv_c1 = - ttt[0][j][i] 
            * (u_wface[0][j  ][i+1] - u_wface[0][j  ][i  ])
                * tarea_r[iblock][j][i] * 2.0; 
        const double adv_c2 = - ttt[0][j][i] 
            * (v_sface[0][j  ][i  ] - v_sface[0][j-1][i  ])
                * tarea_r[iblock][j][i] * 2.0; 

        const double adv_x0 = 
            (ttt[0][j  ][i+1] + ttt[0][j  ][i  ]) * u_wface[0][j  ][i+1] 
                * tarea_r[iblock][j][i]
          - (ttt[0][j  ][i  ] + ttt[0][j  ][i-1]) * u_wface[0][j  ][i  ] 
                * tarea_r[iblock][j][i];
                 
        const double adv_y0 = (
            (ttt[0][j+1][i  ] + ttt[0][j  ][i  ]) * v_sface[0][j  ][i  ]
          - (ttt[0][j  ][i  ] + ttt[0][j-1][i  ]) * v_sface[0][j-1][i  ])
                * tarea_r[iblock][j][i];

        const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

        adv_tt[0][j][i] = adv_xx + adv_yy + adv_zz;
 
        ax[iblock][mtracer][0][j][i] = adv_xx;
        ay[iblock][mtracer][0][j][i] = adv_yy;
        az[iblock][mtracer][0][j][i] = adv_zz;
      }
    }

    for (int k = 1; k < KM-1; ++k) {
      for (int j = 1; j < JMT-1; ++j) {
        for (int i = 1; i < IMT-1; ++i) {
          double adv_xy1;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j  ][i+1] > atmax[k  ][j  ][i+1] || 
              at00[k  ][j  ][i+1] < atmin[k  ][j  ][i+1]) {
         
            adv_xy1 = - (ttt[k][j  ][i+1] - ttt[k][j  ][i  ])
                * fabs(u_wface[k][j  ][i+1]) * tarea_r[iblock][j][i];
          } else {
            adv_xy1 = - dts * (ttt[k][j  ][i+1] - ttt[k][j  ][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i+1], 2)
                    / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
          }

          double adv_xy2;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j  ][i-1] > atmax[k  ][j  ][i-1] || 
              at00[k  ][j  ][i-1] < atmin[k  ][j  ][i-1]) {
         
            adv_xy2 =   (ttt[k][j  ][i  ] - ttt[k][j  ][i-1])
                * fabs(u_wface[k][j  ][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy2 =   dts * (ttt[k][j  ][i  ] - ttt[k][j  ][i-1]) * 2.0
                * tarea_r[iblock][j][i] * pow(u_wface[k][j  ][i  ], 2)
                    / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
          }

          double adv_xy3;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j+1][i  ] > atmax[k  ][j+1][i  ] || 
              at00[k  ][j+1][i  ] < atmin[k  ][j+1][i  ]) {
         
            adv_xy3 = - (ttt[k][j+1][i  ] - ttt[k][j  ][i  ])
                * fabs(v_sface[k][j  ][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy3 = - dts * (ttt[k][j+1][i  ] - ttt[k][j  ][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(v_sface[k][j  ][i  ], 2)
                    / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
          }

          double adv_xy4;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k  ][j-1][i  ] > atmax[k  ][j-1][i  ] || 
              at00[k  ][j-1][i  ] < atmin[k  ][j-1][i  ]) {
         
            adv_xy4 =   (ttt[k][j  ][i  ] - ttt[k][j-1][i  ])
                * fabs(v_sface[k][j-1][i  ]) * tarea_r[iblock][j][i];
          } else {
            adv_xy4 =   dts * (ttt[k][j  ][i  ] - ttt[k][j-1][i  ]) * 2.0
                * tarea_r[iblock][j][i] * pow(v_sface[k][j-1][i  ], 2)
                    / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
          }

          double adv_zb1;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k-1][j  ][i  ] > atmax[k-1][j  ][i  ] || 
              at00[k-1][j  ][i  ] < atmin[k-1][j  ][i  ]) {
  
            adv_zb1 = - 0.5 * fabs(www[k  ][j][i]) * odzp[k]
                * (ttt[k-1][j][i] - ttt[k][j][i]);
          } else {
            adv_zb1 = - 0.5 * odzp[k] * pow(www[k  ][j][i], 2)
                * odzt[k  ] * (ttt[k-1][j][i] - ttt[k  ][j][i]) * dts;
          }

          double adv_zb2;
          if (at00[k  ][j  ][i  ] > atmax[k  ][j  ][i  ] || 
              at00[k  ][j  ][i  ] < atmin[k  ][j  ][i  ] || 
              at00[k+1][j  ][i  ] > atmax[k+1][j  ][i  ] || 
              at00[k+1][j  ][i  ] < atmin[k+1][j  ][i  ]) {
  
            adv_zb2 =   0.5 * fabs(www[k+1][j][i]) * odzp[k]
                * (ttt[k  ][j][i] - ttt[k+1][j][i]);
          } else {
            adv_zb2 =   0.5 * odzp[k] * pow(www[k+1][j][i], 2)
                * odzt[k+1] * (ttt[k  ][j][i] - ttt[k+1][j][i]) * dts;
          }

          const double adv_c1 = - ttt[k][j][i] 
              * (u_wface[k][j  ][i+1] - u_wface[k][j  ][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;
                   
          const double adv_c2 = - ttt[k][j][i] 
              * (v_sface[k][j  ][i  ] - v_sface[k][j-1][i  ]) 
                  * tarea_r[iblock][j][i] * 2.0;

          const double adv_za = 
              0.5 * odzp[k] * www[k  ][j][i] * (ttt[k][j][i] + ttt[k-1][j][i])
            - 0.5 * odzp[k] * www[k+1][j][i] * (ttt[k][j][i] + ttt[k+1][j][i]);

          const double adv_zc = - odzp[k] * ttt[k][j][i]
              * (www[k  ][j][i] - www[k+1][j][i]);

          const double adv_x0 = 
              (ttt[k][j  ][i+1] + ttt[k][j  ][i  ]) * u_wface[k][j  ][i+1] 
                  * tarea_r[iblock][j][i]
            - (ttt[k][j  ][i  ] + ttt[k][j  ][i-1]) * u_wface[k][j  ][i  ] 
                  * tarea_r[iblock][j][i];
                   
          const double adv_y0 = (
              (ttt[k][j+1][i  ] + ttt[k][j  ][i  ]) * v_sface[k][j  ][i  ]
            - (ttt[k][j  ][i  ] + ttt[k][j-1][i  ]) * v_sface[k][j-1][i  ])
                  * tarea_r[iblock][j][i];

          const double adv_xx = - (adv_x0 + adv_xy1 +adv_xy2 + adv_c1);
          const double adv_yy = - (adv_y0 + adv_xy3 +adv_xy4 + adv_c2);
          const double adv_zz = - (adv_za + adv_zb1 +adv_zb2 + adv_zc);

          adv_tt[k][j][i] = adv_xx + adv_yy + adv_zz;

          ax[iblock][mtracer][k][j][i] = adv_xx;
          ay[iblock][mtracer][k][j][i] = adv_yy;
          az[iblock][mtracer][k][j][i] = adv_zz;
        }
      }
    }

    // k = KM - 1
    for (int j = 1; j < JMT-1; ++j) {
      for (int i = 1; i < IMT-1; ++i) {
        double adv_xy1;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i+1] > atmax[KM-1][j  ][i+1] || 
            at00[KM-1][j  ][i+1] < atmin[KM-1][j  ][i+1]) {

          adv_xy1 = - (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ])
              * fabs(u_wface[KM-1][j  ][i+1]) * tarea_r[iblock][j][i];
        } else {
          adv_xy1 = - dts * (ttt[KM-1][j  ][i+1] - ttt[KM-1][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[KM-1][j  ][i+1], 2)
                  / (htw[iblock][j  ][i+1] * hun[iblock][j  ][i+1]);
        }

        double adv_xy2;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i-1] > atmax[KM-1][j  ][i-1] || 
            at00[KM-1][j  ][i-1] < atmin[KM-1][j  ][i-1]) {

          adv_xy2 =   (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1])
              * fabs(u_wface[KM-1][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy2 =   dts * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j  ][i-1]) * 2.0
              * tarea_r[iblock][j][i] * pow(u_wface[KM-1][j  ][i  ], 2)
                  / (htw[iblock][j  ][i  ] * hun[iblock][j  ][i  ]);
        }

        double adv_xy3;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j+1][i  ] > atmax[KM-1][j+1][i  ] || 
            at00[KM-1][j+1][i  ] < atmin[KM-1][j+1][i  ]) {

          adv_xy3 = - (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ])
              * fabs(v_sface[KM-1][j  ][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy3 = - dts * (ttt[KM-1][j+1][i  ] - ttt[KM-1][j  ][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[KM-1][j  ][i  ], 2)
                  / (hts[iblock][j  ][i  ] * hue[iblock][j  ][i  ]);
        }

        double adv_xy4;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-1][j-1][i  ] > atmax[KM-1][j-1][i  ] || 
            at00[KM-1][j-1][i  ] < atmin[KM-1][j-1][i  ]) {

          adv_xy4 =   (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ])
              * fabs(v_sface[KM-1][j-1][i  ]) * tarea_r[iblock][j][i];
        } else {
          adv_xy4 =   dts * (ttt[KM-1][j  ][i  ] - ttt[KM-1][j-1][i  ]) * 2.0
              * tarea_r[iblock][j][i] * pow(v_sface[KM-1][j-1][i  ], 2)
                  / (hts[iblock][j-1][i  ] * hue[iblock][j-1][i  ]);
        }
        double adv_zb1;
        if (at00[KM-1][j  ][i  ] > atmax[KM-1][j  ][i  ] || 
            at00[KM-1][j  ][i  ] < atmin[KM-1][j  ][i  ] || 
            at00[KM-2][j  ][i  ] > atmax[KM-2][j  ][i  ] || 
            at00[KM-2][j  ][i  ] < atmin[KM-2][j  ][i  ]) {
 
          adv_zb1 = - 0.5 * fabs(www[KM-1][j][i]) * odzp[KM-1]
              * (ttt[KM-2][j][i] - ttt[KM-1][j][i]);
        } else {
          adv_zb1 = - 0.5 * odzp[KM-1] * pow(www[KM-1][j][i], 2)
              * odzt[KM-1] * (ttt[KM-2][j][i] - ttt[KM-1][j][i]) * dts;
        }

        const double adv_zb2 = 0.0;

        const double adv_c1 = - ttt[KM-1][j][i] * 
            (u_wface[KM-1][j  ][i+1] - u_wface[KM-1][j  ][i  ])
                * tarea_r[iblock][j][i] * 2.0;
                 
        const double adv_c2 = - ttt[KM-1][j][i] * 
            (v_sface[KM-1][j  ][i  ] - v_sface[KM-1][j-1][i  ])
                * tarea_r[iblock][j][i] * 2.0;

        const double adv_za = 0.5 * odzp[KM-1] * www[KM-1][j][i]
            * (ttt[KM-1][j][i] + ttt[KM-2][j][i]);

        const double adv_zc = - odzp[KM-1] * ttt[KM-1][j][i] * www[KM-1][j][i];

        const double adv_x0 = 
            (ttt[KM-1][j  ][i+1] + ttt[KM-1][j  ][i  ]) 
                * u_wface[KM-1][j  ][i+1] * tarea_r[iblock][j][i]
          - (ttt[KM-1][j  ][i  ] + ttt[KM-1][j  ][i-1])
                * u_wface[KM-1][j  ][i  ] * tarea_r[iblock][j][i];
                 
        const double adv_y0 = (
            (ttt[KM-1][j+1][i  ] + ttt[KM-1][j  ][i  ]) 
                * v_sface[KM-1][j  ][i  ]
          - (ttt[KM-1][j  ][i  ] + ttt[KM-1][j-1][i  ])
                * v_sface[KM-1][j-1][i  ]) * tarea_r[iblock][j][i];

        const double adv_xx = - (adv_x0 + adv_xy1 + adv_xy2 + adv_c1);
        const double adv_yy = - (adv_y0 + adv_xy3 + adv_xy4 + adv_c2);
        const double adv_zz = - (adv_za + adv_zb1 + adv_zb2 + adv_zc);

        adv_tt[KM-1][j][i] = adv_xx + adv_yy + adv_zz;

        ax[iblock][mtracer][KM-1][j][i] = adv_xx;
        ay[iblock][mtracer][KM-1][j][i] = adv_yy;
        az[iblock][mtracer][KM-1][j][i] = adv_zz;
      }
    }
  }
  return ;
}
#endif // LICOM_ENABLE_FORTRAN
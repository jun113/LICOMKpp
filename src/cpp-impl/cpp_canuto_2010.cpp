#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN

#include "../head/cpp_constant_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"

#include <cstdio>
#include <cmath>
#include <vector>
#include <limits>

using CppParamMod::KM;

template<typename T>
inline T sign (const T &x, const T &y) {
  return y >= static_cast<T>(0) ? std::abs(x) : -std::abs(x);
}

void inline cpp_Rf_calc (double &Rf_out, const double &Ri_in, const double &Rrho_in, 
    const double &struct_h, const double &struct_m, const double &r_in) {
  Rf_out = Ri_in * struct_h / struct_m * (1.0 - Rrho_in / r_in);
  return ;
}

void cpp_mixed_layer_tke_calc (double &TKE_out, const double &mld_in, 
    const double &Gm_in, const double &s2_in, const double &lev_in, 
        const double &Rf_in, const double &Rf_inf_in) {
  using CppConstantMod::KARMAN;

  const double l0 = 0.17 * mld_in;
  const double B1 = 21.6;
  const double lB = KARMAN * std::abs(lev_in * l0) / (l0 + KARMAN * std::abs(lev_in));
  const double l  = lB * std::pow(std::pow(1.0 - Rf_in / Rf_inf_in, 4.0), 1.0 / 3.0);

  TKE_out = std::pow(B1, 2) * std::pow(Gm_in, -1.5) * std::pow(l, 2) * std::pow(s2_in, 1.5);
  return ;
}

void cpp_thermocline_mixing_coeff_calc (double &Km_out, 
          double &Kh_out,          double &Ks_out,          double &Kd_out,
    const double &mix_eff_m, const double &mix_eff_h, const double &mix_eff_s, 
    const double &mix_eff_d, const double &n2_in,     const double &lat_in) {
  using CppConstantMod::OMEGA;
  using CppConstantMod::DEGTORAD;
  using CppConstantMod::VERY_SMALL;
  const double N0     = 5.24e-3;
  const double F30    = 2.0 * OMEGA * std::sin(30.0 * DEGTORAD);
  const double Tke_n2 = 0.288e-4;

  double f_lat = std::abs (2.0 * OMEGA * std::sin(lat_in * DEGTORAD) + VERY_SMALL);

  double cond1 = std::max(std::sqrt(n2_in) / f_lat, 1.0);

  double L_lat = (f_lat * std::acosh(cond1)) / (F30 * std::acosh(N0 / F30));

  Km_out = L_lat * Tke_n2 * mix_eff_m;
  Kh_out = L_lat * Tke_n2 * mix_eff_h;
  Ks_out = L_lat * Tke_n2 * mix_eff_s;
  Kd_out = L_lat * Tke_n2 * mix_eff_d;
  return ;
}

void shengjin_calc (const double* array, const int &nmax, double& real_x) {

  const double epsln1 = 1.0e-6;

  const double sa = array[0];
  const double sb = array[1];
  const double sc = array[2];
  const double sd = array[3];

  real_x = 1.0;

  if (std::abs(sa) <= epsln1 && std::abs(sb) <= epsln1) {
    if ((-sd / sc) < 0.0) {
      // printf ("No real positive root at position 1\n");
      return;
    } else {
      // printf ("Only one real positive root exists 1!\n");
      real_x = -sd / sc;
      return;
    }
  }

  double delta;
  double x[3];
  if (std::abs(sa) <= epsln1 && std::abs(sb) >= epsln1) {
    delta = std::pow(sc, 2) - 4.0 * sb * sd;
    if (delta < 0.0) {
      // printf ("No real roots in this secondary order equation 2!\n");
      return;
    } else if (delta >= 0.0) {
      x[0] = (-sc + std::sqrt(delta)) / (2.0 * sb);
      x[1] = (-sc - std::sqrt(delta)) / (2.0 * sb);
      if (x[0] < 0.0 && x[1] < 0.0) {
        // std::cout << "No real positive root at position 2" << std::endl;
        return;
      } else {
        // std::cout << "Two real with one positive at least 2!" << std::endl;
        if (x[0] > 0.0 && x[1] > 0.0) {
          real_x = std::min(x[0], x[1]);
        } else if (x[0] >  0.0 && x[1] <= 0.0) {
          real_x = x[0];
        } else if (x[0] <= 0.0 && x[1]  > 0.0) {
          real_x = x[1];
        }
        return;
      }
    }
  }

  const double A = std::pow(sb, 2) - 3.0 * sa * sc;
  const double B = sb * sc - 9.0 * sa * sd;
  const double C = std::pow(sc, 2) - 3.0 * sb * sd;

  if (std::abs(A) <= epsln1 && std::abs(B) <= epsln1) {
    if ((-sc / sb) < 0.0) {
      // std::cout << "No real positive root at position 3" << std::endl;
      return;
    } else {
      // std::cout << "Only one real positive root exists 3!" << std::endl;
      real_x = -sc / sb;
      return;
    }
  }

  delta = std::pow(B, 2) - 4.0 * A * C;

  double y[2];
  if (delta > 0.0) {
    // TODO
    // y[0] = A * sb + 3.0 * sa * (-B + std::sqrt(delta)) / 2.0;
    // y[1] = A * sb + 3.0 * sa * (-B - std::sqrt(delta)) / 2.0;
    y[0] = A * sb + 1.5 * sa * (-B + std::sqrt(delta));
    y[1] = A * sb + 1.5 * sa * (-B - std::sqrt(delta));
    real_x = (-sb - sign(std::pow(std::abs(y[0]), 1.0 / 3.0), y[0]) 
                  - sign(std::pow(std::abs(y[1]), 1.0 / 3.0), y[1])) 
                      / (3.0 * sa);
    if (real_x < 0.0) {
      real_x = 1.0;
      return;
    } else {
      // std::cout << "Only one real positive root exists 4!" << std::endl;
      // real_x = real_x;
      return;
    }
  }

  const double T = (2.0 * A * sb - 3.0 * sa * B) / (2.0 * std::pow(A, 1.5));
  const double the = std::acos(T);
  x[0] = (-sb - 2.0 * std::sqrt(A) * std::cos(the / 3.0)) / (3.0 * sa);
  x[1] = (-sb + std::sqrt(A) * (std::cos(the / 3.0) 
      + std::sqrt(3.0) * std::sin(the / 3.0))) / (3.0 * sa);
  x[2] = (-sb + std::sqrt(A) * (std::cos(the / 3.0) 
      - std::sqrt(3.0) * std::sin(the / 3.0))) / (3.0 * sa);

  real_x = 500.0;
  for (int k = 0; k < 3; ++k) {
    if (x[k] >= 0.0 && x[k] <= 1000.0) {
      if (x[k] < real_x) {
        real_x = x[k];
      }
    }
  }

  real_x = std::numeric_limits<double>::max();
  for (int i = 0; i < 3; ++i) {
    if (x[i] > 0.0) {
      real_x = std::min (real_x, x[i]);
    }
  }

  real_x = std::min (real_x, 2.0e3);
  real_x = std::max (real_x, 1.0);

  return ;
}

void cpp_dyn_time_scale_calc (const double &Ri, const double &Rrho, 
    const double* pi_n, double &Gm, double* const cubic_coef) {

  double A1, A2, A3, A4, A5, A6;
  // double B1, B2, B3, B4, B5, B6;
  double temp1,  temp2,  temp3,  temp4,  temp5;
  double temp6,  temp7,  temp8,  temp9,  temp10;
  // double temp11, temp12, temp13, temp14, temp15;
  double temp11;

  temp1 = pi_n[0] * pi_n[3] * (pi_n[3] - pi_n[0] * Rrho);
  temp2 = pi_n[1] * (15.0 * pi_n[2] + 7.0) * (std::pow(Rrho, 2) + 1.0);
  temp3 = (14.0 * (pi_n[1] - pi_n[2]) - 15.0 * std::pow(pi_n[2], 2)) * Rrho;
  // hwy_need_check
  temp4 = 150.0 * std::pow(1.0 - Rrho, 3); 

  A1 = temp1 * (temp2 + temp3) / temp4;
  // A1 * temp4 / (1 - Rrho)
  // B1 = pi_n[0] * std::pow(pi_n[3], 2) * (temp2 + temp3) / 150.0;  

  temp1 = pi_n[0] * pi_n[3];
  temp2 = pi_n[1] * (210.0 * pi_n[0] - 150.0 * pi_n[2] + 7.0) 
      * (std::pow(Rrho, 2) + 1.0);
  temp3 = 14.0 * (pi_n[1] - pi_n[2]) * (1.0 + 15.0 * pi_n[0] + 15.0 * pi_n[3]);
  temp4 = 150.0 * pi_n[2] * pi_n[2];
  temp5 = (temp3 + temp4) * Rrho;
  temp6 = 210.0 * pi_n[1] * (pi_n[3] - pi_n[0]);

  // hwy_need_check
  temp7 = 9000.0 * std::pow(1.0 - Rrho, 2); 

  A2 = temp1 * (temp2 + temp5 + temp6) / temp7;
  // B2 = temp1 * (temp2 + temp5 + temp6) / 9000.0;  // A2 * temp7 / (1 - Rrho)

  temp1 = pi_n[0];
  temp2 = 5.0 * pi_n[1] * pi_n[3] * (30.0 * pi_n[2] + 17.0);
  temp3 = pi_n[0] * (15.0 * pi_n[2] + 7.0);
  temp4 = std::pow(Rrho, 2) + 1.0;
  temp5 = temp1 * (temp2 + temp3) * temp4;
  temp6 = -(15.0 * pi_n[2] + 7.0) * (std::pow(pi_n[0], 2) - std::pow(pi_n[3], 2));
  temp7 = 10.0 * pi_n[0] * pi_n[2] * pi_n[3] * (15.0 * pi_n[2] + 17.0);
  temp8 = 15.0 * pi_n[1] * (std::pow(pi_n[0], 2) + std::pow(pi_n[3], 2));
  temp9 = 14.0 * pi_n[0] * pi_n[3] * (1.0 - 10.0 * pi_n[1]);
  temp10 = -(temp7 + temp8 + temp9) * Rrho;
  temp11 = 150.0 * std::pow(1.0 - Rrho, 2);  // hwy_need_check

  A3 = (temp5 + temp6 + temp10) / temp11;
  // B3 = (temp5 + temp6 + temp10) / 150.0;  // A3 * temp11 / (1 - Rrho)

  temp1 = 150.0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
  temp2 = -7.0 * pi_n[0] * (1.0 + 30.0 * pi_n[0]);
  temp3 = (temp1 + temp2) * Rrho;
  temp4 = -150.0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
  temp5 = 7.0 * pi_n[3] * (1.0 + 30.0 * pi_n[3]);
  // temp6 = 9000.0 * std::pow(1.0 - Rrho, 2);  // hwy_need_check
  temp6 = 9000.0 * (1.0 - Rrho);  // hwy_need_check

  A4 = (temp3 + temp4 + temp5) / temp6;
  // B4 = (temp1 + temp2) * (-1);  // B4 * temp6 / (1 - Rrho)

  temp1 = -30.0 * (pi_n[0] * pi_n[2] + pi_n[1] * pi_n[3]);
  temp2 = -17.0 * pi_n[0];
  temp3 = (temp1 + temp2) * Rrho;
  temp4 = 30.0 * (pi_n[0] * pi_n[1] + pi_n[2] * pi_n[3]);
  temp5 = 17.0 * pi_n[3];
  temp6 = 30.0 * (1.0 - Rrho);  // hwy_need_check

  A5 = (temp3 + temp4 + temp5) / temp6;
  // B5 = (temp3 + temp4 + temp5);  // A5 * temp6 / (1 - Rrho)

  A6 = -1.0 / 60.0;
  // B6 = A6;

  cubic_coef[0] = A1 * std::pow(Ri, 3) + A2 * std::pow(Ri, 2);
  cubic_coef[1] = A3 * std::pow(Ri, 2) + A4 * Ri;
  cubic_coef[2] = A5 * Ri + A6;
  cubic_coef[3] = 1.0;

  // call newton_iteration(cubic_coef, 4, Gm);
  shengjin_calc(cubic_coef, 4, Gm);
  // call check_single_value_real_r8_1d("A",[B1,B2,B3,B4,B5,B6],6);
  return ;
}


void cpp_prepare_pi (const double &Ri, const double &Rrho, 
    double* const out_pi) {
  const double a     = 10.0;
  const double Ko    = 1.66;
  const double sig_t = 0.72;
  const double zero  = 1.0e-10;

  const double pi0_1 = 1.0 / (std::sqrt(27.0 / 5.0 * std::pow(Ko, 3)) 
      * (1.0 + 1.0 / sig_t));
  const double pi0_2 = 1.0 / 3.0;
  const double pi0_3 = sig_t;
  const double pi0_4 = pi0_1;
  const double pi0_5 = sig_t;

  if (Ri > zero && Rrho > zero) {
    out_pi[0] = pi0_1 / (1.0 + Ri * Rrho / (a + Rrho));
    out_pi[1] = pi0_2 / (1.0 + Ri) * (1.0 + 2.0 * Ri * Rrho 
        / (1.0 + Rrho * Rrho));
    out_pi[2] = pi0_3;
    out_pi[3] = pi0_4 / (1.0 + Ri / (1.0 + a * Rrho));
    out_pi[4] = pi0_5;
  }

  if (Ri > zero && Rrho <= zero) {
    out_pi[0] = pi0_1 / (1.0 + Ri);
    out_pi[1] = pi0_2;
    out_pi[2] = pi0_3;
    out_pi[3] = pi0_4 / (1.0 + Ri);
    out_pi[4] = pi0_5;
  }

  if (Ri <= zero) {
    out_pi[0] = pi0_1;
    out_pi[1] = pi0_2;
    out_pi[2] = pi0_3;
    out_pi[3] = pi0_4;
    out_pi[4] = pi0_5;
  }
  return ;
}

void cpp_struct_function (const double &Ri, const double &Rrho, const double *pi_n, 
    const double &Gm, double &struct_m, double &struct_h, double &struct_s, 
        double &struct_rho, double &r_out, double &r_out_new) {

  // hwy_need_check
  const double x = Ri * Gm * (1.0 / (1.0 - Rrho));
  const double p = pi_n[3] * pi_n[4] - pi_n[3] * pi_n[1] * (1.0 + Rrho);
  const double q = pi_n[0] * pi_n[1] * (1.0 + Rrho) - pi_n[0] * pi_n[2] * Rrho;

  // hwy_need_check
  const double r     = (pi_n[3] / pi_n[0]) / Rrho * (1.0 + q * x) / (1.0 + p * x);  
  // hwy_need_check
  const double r_new = (pi_n[3] / pi_n[0])        * (1.0 + q * x) / (1.0 + p * x);
  r_out = r;
  r_out_new = r_new;

  const double Ah = pi_n[3] / (1.0 + p * x + pi_n[3] * pi_n[1] * x * (1.0 - 1.0 / r));
  const double As = Ah / r_new;

  // Am1 = 4.0 / 5.0 - (pi_n[3] - pi_n[0] + (pi_n[0] - 1.0 / 150.0) * (1.0 - 1.0 / r)) * x * Ah;
  const double Am1 = 0.8 - (pi_n[3] - pi_n[0] + (pi_n[0] - 1.0 / 150.0) * (1.0 - 1.0 / r)) * x * Ah;
  // const double Am2 = 10.0 + (pi_n[3] - pi_n[0] * Rrho) * x + 1.0 / 50.0 * Gm;
  const double Am2 = 10.0 + (pi_n[3] - pi_n[0] * Rrho) * x + 0.02 * Gm;
  const double Am = Am1 / Am2;

  const double LX = (1.0 - 1.0 / r) * x * Ah;

  const double ratio = 2.0 / 3.0 / (1.0 + 2.0 / 15.0 * LX + 1.0 / 10.0 * Am * Gm);

  struct_m = ratio * Am1 / Am2;
  struct_h = ratio * Ah;
  struct_s = ratio * As;
  struct_rho = (struct_h - struct_s * Rrho) / (1.0 - Rrho);
  return ;
}

void cpp_find_mix_layer_depth (double &mld_out, int &mld_lev_out, 
    const double* den, const double &delta, const double* zlev_in, const int &n) {
  using CppConstantMod::VERY_SMALL;

  double den_m;
  double zlev[n];

  for (int i = 0; i < n; ++i) {
    zlev[i] = std::abs(zlev_in[i]);
  }

  int k;
  for (k = 0; k < n; ++k) {
    if (std::abs(den[k] - den[0]) > delta) {
      den_m    = den[0] - sign(delta, den[0] - den[k]);
      mld_out = zlev[k] + (zlev[k-1] - zlev[k]) * (den_m - den[k]) 
          / (den[k-1] - den[k] + VERY_SMALL);
      mld_lev_out = k + 1;
      break;
    }
  }
  if (k == n) {
    mld_out     = zlev[n - 1];
    mld_lev_out = n;
  }
  return ;
}

void inline cpp_mixing_efficiency (double &mix_eff_var, const double &struct_var, 
    const double &Ri, const double &Gm) {
  mix_eff_var = 0.5 * Ri * Gm * struct_var;
  return ;
}

void cpp_canuto_2010_interface (
          double (&Km_out)[KM],    // wk1
          double (&Kh_out)[KM],    // wk2
          double (&Ks_out)[KM],    // wk3
          double (&Kd_out)[KM],    // wk4
          double (&mld_out),       // amld
    const double (&ts_in)[KM],     // wp1
    const double (&ss_in)[KM],     // wp2
    const double (&rho_in)[KM],    // wp3
    const double (&ri_in)[KM-1],   // wp4
    const double (&rrho_in)[KM-1], // wp5
    const double (&n2_in)[KM-1],   // wp7
    const double (&s2_in)[KM-1],   // wp6
    const double (&lat_in),        // ulat[iblock][j][i] / DegToRad
    const double (&lev_in)[KM],    // wp8
    const int    (&num_lev)) {     /* kmt[iblock][j][i] */

  using CppConstantMod::VERY_SMALL;
  using CppPconstMod::dzp;

  int mld_lev;
  cpp_find_mix_layer_depth (mld_out, mld_lev, rho_in, 0.03, lev_in, num_lev);

  double Ri[KM-1];
  double Rrho[KM-1];
  for (int k = 0; k < KM-1; ++k) {
    Ri[k]   = ri_in[k];
    Rrho[k] = rrho_in[k];
  }

  const double Ri_low  = - 1.0e+10;
  const double Ri_high =   1.0e+10;

  for (int k = 0; k < KM - 1; ++k) {
    if (Ri[k] > Ri_high) {
      Ri[k] = Ri_high;
    } else if (Ri[k] < Ri_low) {
      Ri[k] = Ri_low;
    }
  }

  const double Rrho_bound = 1.0e-3;
  for (int k = 0; k < num_lev - 1; ++k) {
    if (std::abs(Rrho[k] - 1.0) < Rrho_bound) {
      if (Rrho[k] >= 1.0) {
        Rrho[k] = 1.001;
      } else if (Rrho[k] <= 1.0) {
        Rrho[k] = 0.999;
      }
    }
  }
  double out_pi[5];
  double cube[4];
  double Gm[KM - 1];
  double struct_m_inf, struct_h_inf, struct_s_inf, struct_rho_inf;
  double R_inf, Rnew_inf;
  double Rf_inf;
  double struct_m, struct_h, struct_s, struct_rho;
  double R, Rnew;
  double mix_eff_m, mix_eff_h, mix_eff_s, mix_eff_rho;
  double Rf;
  double TKE_mld;
  for (int k = 0; k < num_lev - 1; ++k) {
    cpp_prepare_pi (1.0e10, Rrho[k], out_pi);

    cpp_dyn_time_scale_calc (1.0e10, Rrho[k], out_pi, Gm[k], cube);

    cpp_struct_function(1.0e10, Rrho[k], out_pi, Gm[k], 
        struct_m_inf, struct_h_inf, struct_s_inf, struct_rho_inf, R_inf, Rnew_inf);

    cpp_Rf_calc(Rf_inf, 1.0e10, Rrho[k], struct_h_inf, struct_m_inf, Rnew_inf);

    cpp_prepare_pi(Ri[k], Rrho[k], out_pi);

    cpp_dyn_time_scale_calc(Ri[k], Rrho[k], out_pi, Gm[k], cube);

    cpp_struct_function(Ri[k], Rrho[k], out_pi, Gm[k], 
        struct_m, struct_h, struct_s, struct_rho, R, Rnew);
  
    cpp_mixing_efficiency(mix_eff_m,   struct_m,   Ri[k], Gm[k]);
    cpp_mixing_efficiency(mix_eff_h,   struct_h,   Ri[k], Gm[k]);
    cpp_mixing_efficiency(mix_eff_s,   struct_s,   Ri[k], Gm[k]);
    cpp_mixing_efficiency(mix_eff_rho, struct_rho, Ri[k], Gm[k]);

    cpp_Rf_calc(Rf, Ri[k], Rrho[k], struct_h, struct_m, Rnew);

    if (k <= (mld_lev - 1)) {
      cpp_mixed_layer_tke_calc (TKE_mld, mld_out, Gm[k], 
          s2_in[k], lev_in[k], Rf, Rf_inf);
      Km_out[k] = mix_eff_m   * TKE_mld / (n2_in[k] + VERY_SMALL);
      Kh_out[k] = mix_eff_h   * TKE_mld / (n2_in[k] + VERY_SMALL);
      Ks_out[k] = mix_eff_s   * TKE_mld / (n2_in[k] + VERY_SMALL);
      Kd_out[k] = mix_eff_rho * TKE_mld / (n2_in[k] + VERY_SMALL);
    }

    if (k > (mld_lev - 1)) {
      cpp_thermocline_mixing_coeff_calc(Km_out[k], Kh_out[k], Ks_out[k], Kd_out[k], 
          mix_eff_m, mix_eff_h, mix_eff_s, mix_eff_rho, n2_in[k], lat_in);
    }

    if (k < 3) {
        Km_out[k] = std::max(1.0e-3, Km_out[k]);
        Kh_out[k] = std::max(1.0e-3, Kh_out[k]);
        Ks_out[k] = std::max(1.0e-3, Ks_out[k]);
    } else {
        Km_out[k] = std::max(1.0e-4, Km_out[k]);
        Kh_out[k] = std::max(1.0e-5, Kh_out[k]);
        Ks_out[k] = std::max(1.0e-5, Ks_out[k]);
    }

    Km_out[k] = std::min(1.2e-1, Km_out[k]);
    Kh_out[k] = std::min(1.2e-1, Kh_out[k]);
    Ks_out[k] = std::min(1.2e-1, Ks_out[k]);

    if (n2_in[k] < VERY_SMALL) {
      Kh_out[k] = std::min(8.0e-3 * dzp[k] * dzp[k], 8.0);
      Ks_out[k] = std::min(8.0e-3 * dzp[k] * dzp[k], 8.0);
      Km_out[k] = std::min(8.0e-3 * dzp[k] * dzp[k], 8.0);
    }
  }

  return;
}


// void newton_iteration(const double array[], int nmax, double& out0) {
//     const int max_iter = 1000;
//     const double eps = 1.0e-8;

//     double start0;
//     double x, fx, dfx;
//     bool if_converge;

//     best_initial_guess(array, nmax, start0);
//     shengjin_calc(array, nmax, x, ii, jj);
//     out0 = x;

//    while (!if_converge && counter < max_iter) {
//         fx = 0.0;
//         dfx = 0.0;

//         for (int k = 0; k < nmax; ++k) {
//             fx += array[k] * pow(x, nmax - k);
//         }

//         for (int k = 0; k < nmax - 1; ++k) {
//             dfx += static_cast<double>(nmax - k) * array[k] * pow(x, nmax - 1 - k);
//         }

//         x -= fx / dfx;
//         ++counter;
//         if_converge = std::abs(fx) <= eps;
//     }
//     out0 = x;
//     if (out0 < 0) {
//         std::cerr << "Strong Error: Gm < 0" << std::endl;
//         strat0 = start0 + 1000.0e0;
//     }

//     if (!if_converge) {
//         std::cerr << "After all the iteration times, the program did not converge." << std::endl;
//     }
// }
// // code with return and continue not trans to c++


#endif // LICOM_ENABLE_FORTRAN

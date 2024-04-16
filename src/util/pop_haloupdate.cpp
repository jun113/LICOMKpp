#include "../head/def-undef.h"

#if !defined(LICOM_ENABLE_FORTRAN)

#include "../head/cpp_domain.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_extern_functions.h"
#include "../head/fortran_extern_functions.h"

#ifdef LICOM_ENABLE_KOKKOS
#include "../head/kokkos_dyn_mod.h"
#include "../head/kokkos_work_mod.h"
#include "../head/kokkos_tracer_mod.h"
#include "../head/kokkos_tmp_var.h"

#include "pop_haloupdate.hpp"

#include "Kokkos_Core.hpp"
#endif // LICOM_ENABLE_KOKKOS


#ifdef LICOM_ENABLE_TEST_TIME
#define LICOM_ENABLE_TEST_HALO_UPDATE
#undef  LICOM_ENABLE_TEST_HALO_UPDATE
#endif // LICOM_ENABLE_TEST_TIME

#if (defined (LICOM_ENABLE_KOKKOS)) || (defined (LICOM_ENABLE_CUDA)) || (defined (LICOM_ENABLE_HIP))
// Max size: bclinc2: u, v
// 2 * 4 * km * (imt + jmt)
static double halo_buffer[8 * KM * (IMT + JMT)];
#endif // (define (LICOM_ENABLE_KOKKOS)) || (define (LICOM_ENABLE_CUDA)) || (define (LICOM_ENABLE_HIP))


// For heterogeneous system
void pop_haloupdate_barotr_1(const int &len_k, const int &n_layers) {
  /*
  POP Halo Update: work

  1:  d_halo_buff <- d_work
  2:  h_halo_buff <- d_halo_buff
  3:  work        <- h_halo_buff
  4:  POP Halo Update(work)
  5:  h_halo_buff <- work
  6:  d_halo_buff <- h_halo_buff
  7:  d_work      <- d_halo_buff
  */
#ifdef LICOM_ENABLE_TEST_TIME
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_TIME
  using CppWorkMod::work;
#ifdef LICOM_ENABLE_KOKKOS
  using KokkosWorkMod::p_v_work;

  using Kokkos::parallel_for;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride = total_layers * len_k * IMT;

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), total_layers * len_k * (IMT + JMT));

  static auto v_halo_buffer = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));

  static auto dev = Kokkos::DefaultExecutionSpace();

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_start("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  parallel_for("pop_haloupdate_barotr_work_d2h_i", total_layers * IMT, 
      functor_haloupdate_d2h_i(n_layers, v_halo_buffer, *p_v_work));

  parallel_for("pop_haloupdate_barotr_work_d2h_j", total_layers * JMT, 
      functor_haloupdate_d2h_j(n_layers, v_halo_buffer, *p_v_work));

  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer);

  for (int j = 0; j < n_layers; ++j) {
    for (int i = 0; i < IMT; ++i) {
      work[iblock][j + n_layers][i] = h_v_buffer(i + j * IMT);
    }
  }
  for (int j = n_layers; j < total_layers; ++j) {
    for (int i = 0; i < IMT; ++i) {
      work[iblock][JMT - total_layers - n_layers + j][i] = h_v_buffer(i + j * IMT);
    }
  }
  for (int i = 0; i < n_layers; ++i) {
    for (int j = 0; j < JMT; ++j) {
      work[iblock][j][n_layers + i] = h_v_buffer(stride + i * JMT + j);
    }
  }
  for (int i = n_layers; i < total_layers; ++i) {
    for (int j = 0; j < JMT; ++j) {
      work[iblock][j][IMT - total_layers - n_layers + i] = h_v_buffer(stride + i * JMT + j);
    }
  }

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr memcpy");
  my_time.testTime_start("barotr halo");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  // int errorCode;
  // pop_haloupdate_barotr1_(&errorCode);
  // CppPOPHaloMod::pop_halo_update_2dr8(work[iblock],
  //     CppDomain::POP_haloClinic, 
  //     CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
  //     CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  CppPOPHaloMod::pop_halo_update_2dr8(work[iblock],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr halo");
  my_time.testTime_start("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  for (int j = 0; j < n_layers; ++j) {
    for (int i = 0; i < IMT; ++i) {
      h_v_buffer(i + j * IMT) = work[iblock][j][i];
    }
  }
  for (int j = n_layers; j < total_layers; ++j) {
    for (int i = 0; i < IMT; ++i) {
      h_v_buffer(i + j * IMT) = work[iblock][JMT - total_layers + j][i];
    }
  }

  for (int i = 0; i < n_layers; ++i) {
    for (int j = 0; j < JMT; ++j) {
      h_v_buffer(stride + i * JMT + j) = work[iblock][j][i];
    }
  }
  for (int i = n_layers; i < total_layers; ++i) {
    for (int j = 0; j < JMT; ++j) {
      h_v_buffer(stride + i * JMT + j) = work[iblock][j][IMT - total_layers + i];
    }
  }
  Kokkos::deep_copy(dev, v_halo_buffer, h_v_buffer);
  parallel_for("pop_haloupdate_barotr_work_h2d", total_layers * IMT, 
      functor_haloupdate_h2d_i(n_layers, v_halo_buffer, *p_v_work));
  parallel_for("pop_haloupdate_barotr_work_h2d", total_layers * JMT, 
      functor_haloupdate_h2d_j(n_layers, v_halo_buffer, *p_v_work));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#endif // LICOM_ENABLE_KOKKOS
  return ;
}

void pop_haloupdate_barotr_2(const int &len_k, const int &n_layers) {
/*
POP Halo Update: wka[3], wka[4]

1:  d_halo_buff <- d_wka
2:  h_halo_buff <- d_halo_buff
3:  wka        <- h_halo_buff
4:  POP Halo Update(wka)
5:  h_halo_buff <- work
6:  d_halo_buff <- h_halo_buff
7:  d_wka       <- d_halo_buff
*/
  using CppWorkMod::wka;

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#ifdef LICOM_ENABLE_KOKKOS

  using KokkosWorkMod::p_v_wka;
  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride = len_k * total_layers * IMT;

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), len_k * total_layers * (IMT + JMT));

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, len_k * total_layers * (IMT + JMT)));

  static auto dev = Kokkos::DefaultExecutionSpace();

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_start("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  parallel_for("pop_haloupdate_barotr_wka_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(2, n_layers, v_halo_buffer_sub, *p_v_wka));
  parallel_for("pop_haloupdate_barotr_wka_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(2, len_k, n_layers, v_halo_buffer_sub, *p_v_wka));

  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);

  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][2 + k][j + n_layers][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][2 + k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][2 + k][j][i + n_layers] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][2 + k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr memcpy");
  my_time.testTime_start("barotr halo");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  // int errorCode;
  // pop_haloupdate_barotr2_(&errorCode, &errorCode);
  CppPOPHaloMod::pop_halo_update_2dr8(wka[iblock][2],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
  CppPOPHaloMod::pop_halo_update_2dr8(wka[iblock][3],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
  /*
  CppPOPHaloMod::pop_halo_update(wka[iblock][2], 2, JMT, IMT,
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
  */
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr halo");
  my_time.testTime_start("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) = 
            wka[iblock][2 + k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) = 
            wka[iblock][2 + k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][2 + k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][2 + k][j][IMT - total_layers + i];
      }
    }
  }

  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);

  parallel_for("pop_haloupdate_barotr_wka_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(2, n_layers, v_halo_buffer_sub, *p_v_wka));
  parallel_for("pop_haloupdate_barotr_wka_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(2, len_k, n_layers, v_halo_buffer_sub, *p_v_wka));
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("barotr memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#endif // LICOM_ENABLE_KOKKOS
  return ;
}

void pop_haloupdate_bclinc_2(const int &len_k, const int &n_layers) {

  // u, v
  using CppDynMod::u;
  using CppDynMod::v;
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#ifdef LICOM_ENABLE_KOKKOS
  using KokkosDynMod::p_v_u;
  using KokkosDynMod::p_v_v;

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride   = len_k * total_layers * IMT;
  const int len_arr  = len_k * total_layers * (IMT + JMT);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub_u = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));
  static auto v_halo_buffer_sub_v = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(len_arr, len_arr + total_layers * len_k * (IMT + JMT)));

  // TODO
  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, 2 * total_layers * len_k * (IMT + JMT)));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), 2 * total_layers * len_k * (IMT + JMT));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_start("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  parallel_for("pop_haloupdate_bclinc_u_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub_u, *p_v_u));
  parallel_for("pop_haloupdate_bclinc_u_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub_u, *p_v_u));

  parallel_for("pop_haloupdate_bclinc_v_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub_v, *p_v_v));
  parallel_for("pop_haloupdate_bclinc_v_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub_v, *p_v_v));
  // //----------------
  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);
  // //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        u[iblock][k][j + n_layers][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        u[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        u[iblock][k][j][i + n_layers] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        u[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        v[iblock][k][j + n_layers][i] = 
            h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        v[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        v[iblock][k][j][i + n_layers] = 
            h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        v[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  // //----------------
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc memcpy");
  my_time.testTime_start("bclinc halo");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  // int errorCode;
  // pop_haloupdate_bclinc2_(&errorCode);
  CppPOPHaloMod::pop_halo_update_3dr8(u[iblock],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
  CppPOPHaloMod::pop_halo_update_3dr8(v[iblock],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc halo");
  my_time.testTime_start("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  // //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            u[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            u[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            u[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            u[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i) =
            v[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i) =
            v[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j) =
            v[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j) =
            v[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  // //----------------
  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);
  // //----------------
  parallel_for("pop_haloupdate_bclinc_u_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub_u, *p_v_u));
  parallel_for("pop_haloupdate_bclinc_u_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub_u, *p_v_u));

  parallel_for("pop_haloupdate_bclinc_v_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub_v, *p_v_v));
  parallel_for("pop_haloupdate_bclinc_v_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub_v, *p_v_v));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#endif // LICOM_ENABLE_KOKKOS

  return ;
}

void pop_haloupdate_bclinc_3(const int &len_k, const int &n_layers) {

  using CppWorkMod::wka;
#ifdef LICOM_ENABLE_KOKKOS
  using KokkosWorkMod::p_v_wka;

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride = len_k * total_layers * IMT;

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), total_layers * len_k * (IMT + JMT));

  parallel_for("pop_haloupdate_bclinc_wka_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub, *p_v_wka));
  parallel_for("pop_haloupdate_bclinc_wka_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub, *p_v_wka));
  //----------------
  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][k][j + n_layers][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][k][j][i + n_layers] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  //----------------
  // int errorCode;
  // pop_haloupdate_bclinc3_(&errorCode);
  CppPOPHaloMod::pop_halo_update_3dr8(wka[iblock],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_SW_CORNER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_VECTOR);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            wka[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            wka[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  //----------------
  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);
  //----------------
  parallel_for("pop_haloupdate_bclinc_wka_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub, *p_v_wka));

  parallel_for("pop_haloupdate_bclinc_wka_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub, *p_v_wka));

#endif // LICOM_ENABLE_KOKKOS

  return ;
}


/*
void pop_haloupdate_bclinc_33 (const int &len_k, const int &n_layers) {

  using CppWorkMod::wka;
  using CppWorkMod::wkb;
  int errorCode;

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#ifdef LICOM_ENABLE_KOKKOS
  using KokkosWorkMod::p_v_wka;
  using KokkosWorkMod::p_v_wkb;

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride = len_k * total_layers * IMT;
  const int len_arr  = len_k * total_layers * (IMT + JMT);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, 2 * total_layers * len_k * (IMT + JMT)));
  static auto v_halo_buffer_sub_wka = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));
  static auto v_halo_buffer_sub_wkb = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(len_arr, len_arr + total_layers * len_k * (IMT + JMT)));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), 2 * total_layers * len_k * (IMT + JMT));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_start("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  parallel_for("pop_haloupdate_bclinc_wka_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub_wka, *p_v_wka));
  parallel_for("pop_haloupdate_bclinc_wka_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub_wka, *p_v_wka));
  parallel_for("pop_haloupdate_bclinc_wkb_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub_wkb, *p_v_wkb));
  parallel_for("pop_haloupdate_bclinc_wkb_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub_wkb, *p_v_wkb));
  //----------------
  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][k][j + n_layers][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wka[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][k][j][i + n_layers] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wka[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wkb[iblock][k][j + n_layers][i] = 
            h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        wkb[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wkb[iblock][k][j][i + n_layers] = 
            h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        wkb[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  //----------------
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc memcpy");
  my_time.testTime_start("bclinc halo");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  pop_haloupdate_bclinc33_(&errorCode);
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc halo");
  my_time.testTime_start("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            wka[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            wka[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            wka[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i) =
            wkb[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(len_arr + k * total_layers * IMT + j * IMT + i) =
            wkb[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j) =
            wkb[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(len_arr + stride + k * total_layers * JMT + i * JMT + j) =
            wkb[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  //----------------
  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);
  //----------------
  parallel_for("pop_haloupdate_bclinc_wka_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub_wka, *p_v_wka));

  parallel_for("pop_haloupdate_bclinc_wka_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub_wka, *p_v_wka));
  parallel_for("pop_haloupdate_bclinc_wkb_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub_wkb, *p_v_wkb));

  parallel_for("pop_haloupdate_bclinc_wkb_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub_wkb, *p_v_wkb));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("bclinc memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#endif // LICOM_ENABLE_KOKKOS

  return ;
}
*/

void pop_haloupdate_tracer_1(const int &len_k, const int &n_layers) {

  // net
  using CppTracerMod::net;
#ifdef LICOM_ENABLE_KOKKOS

  using KokkosTracerMod::p_v_net;
  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int total_layers = n_layers << 1;
  const int stride = len_k * total_layers * IMT;

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), total_layers * len_k * (IMT + JMT));

  parallel_for("pop_haloupdate_tracer_net_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_d2h_i(n_layers, v_halo_buffer_sub, *p_v_net));
  parallel_for("pop_haloupdate_tracer_net_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_d2h_j(len_k, n_layers, v_halo_buffer_sub, *p_v_net));
  //----------------
  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        net[iblock][k][j + n_layers][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        net[iblock][k][JMT - total_layers - n_layers + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        net[iblock][k][j][i + n_layers] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        net[iblock][k][j][IMT - total_layers - n_layers + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
  //----------------
  int errorCode;
  pop_haloupdate_tracer1_(&errorCode);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            net[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            net[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            net[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            net[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  //----------------
  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);
  //----------------
  parallel_for("pop_haloupdate_tracer_net_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          functor_haloupdate_h2d_i(n_layers, v_halo_buffer_sub, *p_v_net));

  parallel_for("pop_haloupdate_tracer_net_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          functor_haloupdate_h2d_j(len_k, n_layers, v_halo_buffer_sub, *p_v_net));

#endif // LICOM_ENABLE_KOKKOS

  return ;
}

void pop_haloupdate_tracer_2(const int &len_k) {


#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  using MyTest::my_time;
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  using CppDynMod::vtl;
#ifdef LICOM_ENABLE_KOKKOS

  using KokkosDynMod::p_v_vtl;
  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  const int iblock = 0;
  const int n_layers = 3;
  const int total_layers = 6;
  const int stride = len_k * total_layers * IMT;

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, total_layers * len_k * (IMT + JMT)));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(&(halo_buffer[0]), total_layers * len_k * (IMT + JMT));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_start("tracer memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

  parallel_for("pop_haloupdate_tracer_vtl_d2h_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          FuncHaloUpdateD2HTracVtlI(n_layers, v_halo_buffer_sub, *p_v_vtl));
  parallel_for("pop_haloupdate_tracer_vtl_d2h_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          FuncHaloUpdateD2HTracVtlJ(len_k, n_layers, v_halo_buffer_sub, *p_v_vtl));
  //----------------
  Kokkos::deep_copy(dev, h_v_buffer, v_halo_buffer_sub);
  //----------------
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        vtl[iblock][k][j + 2][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        vtl[iblock][k][JMT - total_layers - 2 + j][i] = 
            h_v_buffer(k * total_layers * IMT + j * IMT + i);
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        vtl[iblock][k][j][i + 2] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        vtl[iblock][k][j][IMT - total_layers - 2 + i] = 
            h_v_buffer(stride + k * total_layers * JMT + i * JMT + j);
      }
    }
  }
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("tracer memcpy");
  my_time.testTime_start("tracer halo");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  //----------------
  // int errorCode;
  // pop_haloupdate_tracer2_(&errorCode);
  CppPOPHaloMod::pop_halo_update_3dr8(vtl[iblock],
      CppDomain::POP_haloClinic_C, 
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  //----------------
#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("tracer halo");
  my_time.testTime_start("tracer memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE
  for (int k = 0; k < len_k; ++k) {
    for (int j = 0; j < n_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            vtl[iblock][k][j][i];
      }
    }
    for (int j = n_layers; j < total_layers; ++j) {
      for (int i = 0; i < IMT; ++i) {
        h_v_buffer(k * total_layers * IMT + j * IMT + i) =
            vtl[iblock][k][JMT - total_layers + j][i];
      }
    }
  }
  for (int k = 0; k < len_k; ++k) {
    for (int i = 0; i < n_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            vtl[iblock][k][j][i];
      }
    }
    for (int i = n_layers; i < total_layers; ++i) {
      for (int j = 0; j < JMT; ++j) {
        h_v_buffer(stride + k * total_layers * JMT + i * JMT + j) =
            vtl[iblock][k][j][IMT - total_layers + i];
      }
    }
  }
  //----------------
  Kokkos::deep_copy(dev, v_halo_buffer_sub, h_v_buffer);
  //----------------
  parallel_for("pop_haloupdate_tracer_vtl_h2d_i", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * IMT, len_k}), 
          FuncHaloUpdateH2DTracVtlI(n_layers, v_halo_buffer_sub, *p_v_vtl));

  parallel_for("pop_haloupdate_tracer_vtl_h2d_j", 
      MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {total_layers * JMT, len_k}), 
          FuncHaloUpdateH2DTracVtlJ(len_k, n_layers, v_halo_buffer_sub, *p_v_vtl));

#ifdef LICOM_ENABLE_TEST_HALO_UPDATE
  my_time.testTime_stop("tracer memcpy");
#endif // LICOM_ENABLE_TEST_HALO_UPDATE

#endif // LICOM_ENABLE_KOKKOS
}

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
// void gpu_get_halo_transpose_double (const ViewDouble4D &viewSrc, double* const arrObj,
//     const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

//   using Kokkos::parallel_for;
//   using Kokkos::MDRangePolicy;

//   int start_b[4], end_b[4];
//   int start_c[4], end_c[4];

//   start_b[0] = startLayer;
//   start_b[1] = startLayer + lenLayer;
//   start_b[2] = startLayer + lenLayer;
//   start_b[3] = lenB - startLayer - lenLayer;
//   end_b[0]   = startLayer + lenLayer;
//   end_b[1]   = lenB - startLayer - lenLayer;
//   end_b[2]   = lenB - startLayer - lenLayer;
//   end_b[3]   = lenB - startLayer;
//   start_c[0] = startLayer;
//   start_c[1] = startLayer;
//   start_c[2] = lenC - startLayer - lenLayer;
//   start_c[3] = startLayer;
//   end_c[0]   = lenC - startLayer;
//   end_c[1]   = startLayer + lenLayer;
//   end_c[2]   = lenC - startLayer;
//   end_c[3]   = lenC - startLayer;

//   int num_block[4];
//   num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
//   num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);
//   num_block[2] = (end_b[2] - start_b[2]) * (end_c[2] - start_c[2]);
//   num_block[3] = (end_b[3] - start_b[3]) * (end_c[3] - start_c[3]);

//   const int total_num_block = num_block[0] + num_block[1]
//                             + num_block[2] + num_block[3];

//   const int max_num_block = std::max(num_block[0], num_block[1]);

//   static auto dev = Kokkos::DefaultExecutionSpace();

//   static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
//       Kokkos::make_pair(0, lenA * total_num_block));

//   static Kokkos::View<double *, Layout, 
//       Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
//           h_v_buffer(halo_buffer, lenA * total_num_block);

//   const int team_size   = 128;
//   const int league_size = (lenA * max_num_block + team_size - 1) / team_size;
//   // const int league_size = lenA * ((max_num_block + team_size - 1) / team_size);
 
//   using ScratchView = Kokkos::View<double*, Kokkos::DefaultExecutionSpace::scratch_memory_space>;
//   // size_t shmem_size = ScratchView::shmem_size (Kokkos::AUTO());
 
//   parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
//       league_size, team_size),
//           FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
//               startLayer, lenLayer, lenA, lenB, lenC, max_num_block, CppParamMod::mytid));

//   Kokkos::deep_copy (dev, h_v_buffer, v_halo_buffer_sub);

//   const int stride_obj = lenC * lenA;

//   int len_block_b = end_b[0] - start_b[0];
//   int len_block_c = end_c[0] - start_c[0];
//   for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
//     for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
//       for (int idx_a = 0; idx_a < lenA; ++idx_a) {
//         arrObj[(start_b[0]+idx_b)*stride_obj + (start_c[0]+idx_c)*lenA + idx_a] =
//             h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c);
//       }
//     }
//   }

//   int start_buff = lenA * num_block[0];
//   len_block_b = end_b[1] - start_b[1];
//   len_block_c = end_c[1] - start_c[1];
//   for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
//     for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
//       for (int idx_a = 0; idx_a < lenA; ++idx_a) {
//         arrObj[(start_b[1]+idx_b)*stride_obj + (start_c[1]+idx_c)*lenA + idx_a] =
//             h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
//       }
//     }
//   }

//   start_buff += (lenA * num_block[1]);
//   len_block_b = end_b[2] - start_b[2];
//   len_block_c = end_c[2] - start_c[2];
//   for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
//     for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
//       for (int idx_a = 0; idx_a < lenA; ++idx_a) {
//         arrObj[(start_b[2]+idx_b)*stride_obj + (start_c[2]+idx_c)*lenA + idx_a] =
//             h_v_buffer(start_buff + idx_a*num_block[2] + idx_b*len_block_c + idx_c);
//       }
//     }
//   }

//   start_buff += (lenA * num_block[2]);
//   len_block_b = end_b[3] - start_b[3];
//   len_block_c = end_c[3] - start_c[3];
//   for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
//     for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
//       for (int idx_a = 0; idx_a < lenA; ++idx_a) {
//         arrObj[(start_b[3]+idx_b)*stride_obj + (start_c[3]+idx_c)*lenA + idx_a] =
//             h_v_buffer(start_buff + idx_a*num_block[3] + idx_b*len_block_c + idx_c);
//       }
//     }
//   }

//   return ;
// }

void gpu_get_halo_transpose_double (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] << 1) 
                            + (num_block[1] << 1);

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  // const int team_size   = 256;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;
  // const int league_size = lenA * ((max_num_block + team_size - 1) / team_size);
 
  // using ScratchView = Kokkos::View<double*, Kokkos::DefaultExecutionSpace::scratch_memory_space>;
  // size_t shmem_size = ScratchView::shmem_size (1056);

  parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));
 
  // parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
  //     league_size, team_size).set_scratch_size(0,Kokkos::PerTeam(shmem_size)),
  //         FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
  //             startLayer, lenLayer, lenA, lenB, lenC, max_num_block));

  // parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
  //     league_size, team_size).set_scratch_size(0,Kokkos::PerTeam(shmem_size)),
  //         FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
  //             startLayer, lenLayer, lenA, lenB, lenC));

  Kokkos::deep_copy (dev, h_v_buffer, v_halo_buffer_sub);

  const int stride_obj = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  // int stride_buff = len_block_c * lenA;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[0]+idx_b)*stride_obj + (start_c[0]+idx_c)*lenA + idx_a] =
            h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c);
        // arrObj[(start_b[0]+idx_b)*stride_obj + (start_c[0]+idx_c)*lenA + idx_a] =
        //     h_v_buffer(idx_b*stride_buff + idx_c*lenA + idx_a);
      }
    }
  }

  int start_buff = lenA * num_block[0];
  // int start_buff = lenA * len_block_c;
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  // stride_buff = len_block_c * lenA;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[1]+idx_b)*stride_obj + (start_c[1]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
        // arrObj[(start_b[1]+idx_b)*stride_obj + (start_c[1]+idx_c)*lenA + idx_a] =
        //     h_v_buffer(start_buff + idx_b*stride_buff + idx_c*lenA + idx_a);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  // start_buff += (lenA * len_block_c);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  // stride_buff = len_block_c * lenA;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[2]+idx_b)*stride_obj + (start_c[2]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
        // arrObj[(start_b[2]+idx_b)*stride_obj + (start_c[2]+idx_c)*lenA + idx_a] =
        //     h_v_buffer(start_buff + idx_b*stride_buff + idx_c*lenA + idx_a);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  // start_buff += (lenA * len_block_c);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  // stride_buff = len_block_c * lenA;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[3]+idx_b)*stride_obj + (start_c[3]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c);
        // arrObj[(start_b[3]+idx_b)*stride_obj + (start_c[3]+idx_c)*lenA + idx_a] =
        //     h_v_buffer(start_buff + idx_b*stride_buff + idx_c*lenA + idx_a);
      }
    }
  }

  return ;
}

void gpu_put_halo_transpose_double (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] + num_block[1]) << 1;

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;
  // const int league_size = lenA * ((max_num_block + team_size - 1) / team_size);
  // const int league_size = (lenA + team_size - 1) / team_size;
 
  // using ScratchView = Kokkos::View<double*, Kokkos::DefaultExecutionSpace::scratch_memory_space>;
  // size_t shmem_size = ScratchView::shmem_size (team_size);

  const int stride_src = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[0]+idx_b)*stride_src + (start_c[0]+idx_c)*lenA + idx_a];
      }
    }
  }

  int start_buff = lenA * num_block[0];
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[1]+idx_b)*stride_src + (start_c[1]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[2]+idx_b)*stride_src + (start_c[2]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[3]+idx_b)*stride_src + (start_c[3]+idx_c)*lenA + idx_a];
      }
    }
  }

  Kokkos::deep_copy (dev, v_halo_buffer_sub, h_v_buffer);

  // parallel_for ("h2d_pop_halo_trans_k", Kokkos::TeamPolicy<>(
  //     league_size, team_size).set_scratch_size(0, Kokkos::PerTeam(shmem_size)),
  //         FuncPutHaloTransDouble(v_halo_buffer_sub, viewDst,
  //             startLayer, lenLayer, lenA, lenB, lenC, max_num_block));
  parallel_for ("h2d_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncPutHaloTransDouble(v_halo_buffer_sub, viewDst,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));
  return ;
}
void gpu_get_halo_transpose_bclinc (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] << 1) 
                            + (num_block[1] << 1);

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;

  parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));

  Kokkos::deep_copy (dev, h_v_buffer, v_halo_buffer_sub);

  const int stride_obj = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[0]+idx_b)*stride_obj + (start_c[0]+idx_c)*lenA + idx_a] =
            h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c);
      }
    }
  }

  int start_buff = lenA * num_block[0];
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[1]+idx_b)*stride_obj + (start_c[1]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[2]+idx_b)*stride_obj + (start_c[2]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[3]+idx_b)*stride_obj + (start_c[3]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c);
      }
    }
  }

  return ;
}

void gpu_put_halo_transpose_bclinc (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] + num_block[1]) << 1;

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;

  const int stride_src = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[0]+idx_b)*stride_src + (start_c[0]+idx_c)*lenA + idx_a];
      }
    }
  }

  int start_buff = lenA * num_block[0];
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[1]+idx_b)*stride_src + (start_c[1]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[2]+idx_b)*stride_src + (start_c[2]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[3]+idx_b)*stride_src + (start_c[3]+idx_c)*lenA + idx_a];
      }
    }
  }

  Kokkos::deep_copy (dev, v_halo_buffer_sub, h_v_buffer);

  parallel_for ("h2d_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncPutHaloTransDouble(v_halo_buffer_sub, viewDst,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));
  return ;
}
void gpu_get_halo_transpose_tracer (const ViewDouble4D &viewSrc, double* const arrObj,
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] << 1) 
                            + (num_block[1] << 1);

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;

  parallel_for ("d2h_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncGetHaloTransDouble(viewSrc, v_halo_buffer_sub,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));

  Kokkos::deep_copy (dev, h_v_buffer, v_halo_buffer_sub);

  const int stride_obj = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[0]+idx_b)*stride_obj + (start_c[0]+idx_c)*lenA + idx_a] =
            h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c);
      }
    }
  }

  int start_buff = lenA * num_block[0];
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[1]+idx_b)*stride_obj + (start_c[1]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[2]+idx_b)*stride_obj + (start_c[2]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c);
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        arrObj[(start_b[3]+idx_b)*stride_obj + (start_c[3]+idx_c)*lenA + idx_a] =
            h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c);
      }
    }
  }

  return ;
}

void gpu_put_halo_transpose_tracer (double* const arrSrc, const ViewDouble4D &viewDst, 
    const int &startLayer, const int &lenLayer, const int &lenA, const int &lenB, const int &lenC) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  int start_b[4], end_b[4];
  int start_c[4], end_c[4];

  start_b[0] = startLayer;
  start_b[1] = startLayer + lenLayer;
  start_b[2] = startLayer + lenLayer;
  start_b[3] = lenB - startLayer - lenLayer;
  end_b[0]   = startLayer + lenLayer;
  end_b[1]   = lenB - startLayer - lenLayer;
  end_b[2]   = lenB - startLayer - lenLayer;
  end_b[3]   = lenB - startLayer;
  start_c[0] = startLayer;
  start_c[1] = startLayer;
  start_c[2] = lenC - startLayer - lenLayer;
  start_c[3] = startLayer;
  end_c[0]   = lenC - startLayer;
  end_c[1]   = startLayer + lenLayer;
  end_c[2]   = lenC - startLayer;
  end_c[3]   = lenC - startLayer;

  int num_block[2];
  num_block[0] = (end_b[0] - start_b[0]) * (end_c[0] - start_c[0]);
  num_block[1] = (end_b[1] - start_b[1]) * (end_c[1] - start_c[1]);

  const int total_num_block = (num_block[0] + num_block[1]) << 1;

  const int max_num_block = std::max(num_block[0], num_block[1]);

  static auto dev = Kokkos::DefaultExecutionSpace();

  static auto v_halo_buffer_sub = Kokkos::subview(*KokkosTmpVar::p_v_halo_buffer,
      Kokkos::make_pair(0, lenA * total_num_block));

  static Kokkos::View<double *, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_buffer(halo_buffer, lenA * total_num_block);

  const int team_size   = 128;
  const int league_size = (lenA * max_num_block + team_size - 1) / team_size;

  const int stride_src = lenC * lenA;

  int len_block_b = lenLayer;
  int len_block_c = end_c[0] - start_c[0];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[0]+idx_b)*stride_src + (start_c[0]+idx_c)*lenA + idx_a];
      }
    }
  }

  int start_buff = lenA * num_block[0];
  len_block_b = end_b[1] - start_b[1];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[1]+idx_b)*stride_src + (start_c[1]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = end_b[2] - start_b[2];
  len_block_c = lenLayer;
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[1] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[2]+idx_b)*stride_src + (start_c[2]+idx_c)*lenA + idx_a];
      }
    }
  }

  start_buff += (lenA * num_block[1]);
  len_block_b = lenLayer;
  len_block_c = end_c[3] - start_c[3];
  for (int idx_b = 0; idx_b < len_block_b; ++idx_b) {
    for (int idx_c = 0; idx_c < len_block_c; ++idx_c) {
      for (int idx_a = 0; idx_a < lenA; ++idx_a) {
        h_v_buffer(start_buff + idx_a*num_block[0] + idx_b*len_block_c + idx_c) = 
            arrSrc[(start_b[3]+idx_b)*stride_src + (start_c[3]+idx_c)*lenA + idx_a];
      }
    }
  }

  Kokkos::deep_copy (dev, v_halo_buffer_sub, h_v_buffer);

  parallel_for ("h2d_pop_halo_trans_k", Kokkos::TeamPolicy<>(
      league_size, team_size),
          FuncPutHaloTransDouble(v_halo_buffer_sub, viewDst,
              startLayer, lenLayer, lenA, lenB, lenC, max_num_block));
  return ;
}
#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

#endif // !defined(LICOM_ENABLE_FORTRAN)

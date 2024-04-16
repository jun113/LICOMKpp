#include "../head/def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include "kokkos_jra_daily.hpp"

#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/kokkos_forc_mod.h"
#include "../head/kokkos_pconst_mod.h"

#include "../head/cpp_extern_functions.h"
//=========
#include "../head/cpp_forc_mod.h"
#include "../head/cpp_param_mod.h"
#include "../head/cpp_pconst_mod.h"
#include "../head/cpp_pop_halo_mod.hpp"
#include "../head/cpp_pop_grid_horz_mod.h"
#include "../head/cpp_msg_mod.h"

void kokkos_init_jra_daily_low () {
  using KokkosForcMod::p_v_t10;
  using KokkosForcMod::p_v_u10;
  using KokkosForcMod::p_v_v10;
  using KokkosForcMod::p_v_slp;
  using KokkosForcMod::p_v_q10;
  using KokkosForcMod::p_v_swhf;
  using KokkosForcMod::p_v_lwhf;
  using KokkosForcMod::p_v_precr;
  using KokkosForcMod::p_v_precs;
  using KokkosForcMod::p_v_rf;
  using KokkosForcMod::p_v_si;
  using KokkosForcMod::p_v_buffer;
  using KokkosPconstMod::p_v_s_lon;
  using KokkosPconstMod::p_v_s_lat;

  CppForcMod::t10   = new float[S_IMT * S_JMT];
  CppForcMod::u10   = new float[S_IMT * S_JMT];
  CppForcMod::v10   = new float[S_IMT * S_JMT];
  CppForcMod::slp   = new float[S_IMT * S_JMT];
  CppForcMod::q10   = new float[S_IMT * S_JMT];
  CppForcMod::swhf  = new float[S_IMT * S_JMT];
  CppForcMod::lwhf  = new float[S_IMT * S_JMT];
  CppForcMod::precr = new float[S_IMT * S_JMT];
  CppForcMod::precs = new float[S_IMT * S_JMT];
  CppForcMod::rf    = new float[S_IMT * S_JMT];
  CppForcMod::si    = new float[S_IMT * S_JMT];

  p_v_t10   = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_u10   = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_v10   = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_slp   = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_q10   = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_swhf  = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_lwhf  = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_precr = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_precs = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_rf    = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
  p_v_si    = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));

  // jra daily
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  CppForcMod::tsa3    = new double [JMT * IMT];
  CppForcMod::wspd3   = new double [JMT * IMT];
  CppForcMod::wspdu3  = new double [JMT * IMT];
  CppForcMod::wspdv3  = new double [JMT * IMT];
  CppForcMod::psa3    = new double [JMT * IMT];
  CppForcMod::qar3    = new double [JMT * IMT];
  CppForcMod::swv3    = new double [JMT * IMT];
  CppForcMod::lwv3    = new double [JMT * IMT];
  CppForcMod::rain3   = new double [JMT * IMT];
  CppForcMod::snow3   = new double [JMT * IMT];
  CppForcMod::runoff3 = new double [JMT * IMT];
  CppForcMod::seaice3 = new double [JMT * IMT];
  new (p_v_t10)   ViewFloat3D("pointer_view_t10",   1, S_JMT, S_IMT);
  new (p_v_u10)   ViewFloat3D("pointer_view_u10",   1, S_JMT, S_IMT);
  new (p_v_v10)   ViewFloat3D("pointer_view_v10",   1, S_JMT, S_IMT);
  new (p_v_slp)   ViewFloat3D("pointer_view_slp",   1, S_JMT, S_IMT);
  new (p_v_q10)   ViewFloat3D("pointer_view_q10",   1, S_JMT, S_IMT);
  new (p_v_swhf)  ViewFloat3D("pointer_view_swhf",  1, S_JMT, S_IMT);
  new (p_v_lwhf)  ViewFloat3D("pointer_view_lwhf",  1, S_JMT, S_IMT);
  new (p_v_precr) ViewFloat3D("pointer_view_precr", 1, S_JMT, S_IMT);
  new (p_v_precs) ViewFloat3D("pointer_view_precs", 1, S_JMT, S_IMT);
  new (p_v_rf)    ViewFloat3D("pointer_view_rf",    1, S_JMT, S_IMT);
  new (p_v_si)    ViewFloat3D("pointer_view_si",    1, S_JMT, S_IMT);
#else
  new (p_v_t10)   ViewFloat3D(CppForcMod::t10,   1, S_JMT, S_IMT);
  new (p_v_u10)   ViewFloat3D(CppForcMod::u10,   1, S_JMT, S_IMT);
  new (p_v_v10)   ViewFloat3D(CppForcMod::v10,   1, S_JMT, S_IMT);
  new (p_v_slp)   ViewFloat3D(CppForcMod::slp,   1, S_JMT, S_IMT);
  new (p_v_q10)   ViewFloat3D(CppForcMod::q10,   1, S_JMT, S_IMT);
  new (p_v_swhf)  ViewFloat3D(CppForcMod::swhf,  1, S_JMT, S_IMT);
  new (p_v_lwhf)  ViewFloat3D(CppForcMod::lwhf,  1, S_JMT, S_IMT);
  new (p_v_precr) ViewFloat3D(CppForcMod::precr, 1, S_JMT, S_IMT);
  new (p_v_precs) ViewFloat3D(CppForcMod::precs, 1, S_JMT, S_IMT);
  new (p_v_rf)    ViewFloat3D(CppForcMod::rf,    1, S_JMT, S_IMT);
  new (p_v_si)    ViewFloat3D(CppForcMod::si,    1, S_JMT, S_IMT);
#endif

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  ViewDouble2D::HostMirror h_v_s_lon = create_mirror_view(*p_v_s_lon);
  ViewDouble2D::HostMirror h_v_s_lat = create_mirror_view(*p_v_s_lat);

  ViewDouble1D::HostMirror h_v_buffer = create_mirror_view(*p_v_buffer);

  // update s_lon s_lat
  read_jra("runoff_all.2016_daily.nc", h_v_buffer.data(),
      h_v_s_lon.data(), h_v_s_lat.data());
  Kokkos::deep_copy(*p_v_s_lon, h_v_s_lon);
  Kokkos::deep_copy(*p_v_s_lat, h_v_s_lat);

#else
  read_jra("runoff_all.2016_daily.nc", CppForcMod::buffer,
      (*p_v_s_lon).data(), (*p_v_s_lat).data());
#endif

  return ;
}

void kokkos_init_jra_daily_high() {

  using KokkosForcMod::p_v_t10;
  using KokkosForcMod::p_v_u10;
  using KokkosForcMod::p_v_v10;
  using KokkosForcMod::p_v_slp;
  using KokkosForcMod::p_v_q10;
  using KokkosForcMod::p_v_swhf;
  using KokkosForcMod::p_v_lwhf;
  using KokkosForcMod::p_v_precr;
  using KokkosForcMod::p_v_precs;
  using KokkosForcMod::p_v_rf;
  using KokkosForcMod::p_v_si;
  using KokkosForcMod::p_v_buffer;
  using KokkosPconstMod::p_v_s_lon;
  using KokkosPconstMod::p_v_s_lat;

  if (CppMsgMod::nproc < 11) {
    if (CppParamMod::mytid == 0) {
    printf("The number of processes must be greater than 11 at high-resolution, currently %d !\n", 
        CppMsgMod::nproc);
  }
  exit(EXIT_FAILURE);
  }
  const int irec = CppPconstMod::number_day;
  if (CppParamMod::mytid == 0) {
  printf("jra daily start irec = %d\n", irec);
  }

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  CppForcMod::tsa3    = new double [JMT * IMT];
  CppForcMod::wspd3   = new double [JMT * IMT];
  CppForcMod::wspdu3  = new double [JMT * IMT];
  CppForcMod::wspdv3  = new double [JMT * IMT];
  CppForcMod::psa3    = new double [JMT * IMT];
  CppForcMod::qar3    = new double [JMT * IMT];
  CppForcMod::swv3    = new double [JMT * IMT];
  CppForcMod::lwv3    = new double [JMT * IMT];
  CppForcMod::rain3   = new double [JMT * IMT];
  CppForcMod::snow3   = new double [JMT * IMT];
  CppForcMod::runoff3 = new double [JMT * IMT];
  CppForcMod::seaice3 = new double [JMT * IMT];
  ViewDouble2D::HostMirror h_v_s_lon = create_mirror_view(*p_v_s_lon);
  ViewDouble2D::HostMirror h_v_s_lat = create_mirror_view(*p_v_s_lat);
#endif

  // ! read in jra data
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  if (CppParamMod::mytid == 0) {
    p_v_rf = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
    new (p_v_rf) ViewFloat3D("pointer_view_rf", 365, S_JMT, S_IMT);
    ViewFloat3D::HostMirror h_v_rf = create_mirror_view(*p_v_rf);
    ViewDouble1D::HostMirror h_v_buffer = create_mirror_view(*p_v_buffer);
    // update s_lon s_lat
    read_jra(irec, 365, "runoff_all.2016_daily.nc", h_v_buffer.data(), 
        h_v_s_lon.data(), h_v_s_lat.data(), h_v_rf.data());
    Kokkos::deep_copy(*p_v_rf, h_v_rf);
  } else if (CppParamMod::mytid == 1) {
      p_v_t10 = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_t10) ViewFloat3D("pointer_view_t10", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_t10 = create_mirror_view(*p_v_t10);
    read_jra(irec, 365, "t_10.2016_daily.nc", h_v_t10.data());
    Kokkos::deep_copy(*p_v_t10, h_v_t10);
  } else if (CppParamMod::mytid == 2) {
      p_v_u10 = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_u10) ViewFloat3D("pointer_view_u10", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_u10 = create_mirror_view(*p_v_u10);
    read_jra(irec, 365, "u_10.2016_daily.nc", h_v_u10.data());
    Kokkos::deep_copy(*p_v_u10, h_v_u10);
  } else if (CppParamMod::mytid == 3) {
      p_v_v10 = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_v10) ViewFloat3D("pointer_view_v10", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_v10 = create_mirror_view(*p_v_v10);
    read_jra(irec, 365, "v_10.2016_daily.nc", h_v_v10.data());
    Kokkos::deep_copy(*p_v_v10, h_v_v10);
  } else if (CppParamMod::mytid == 4) {
      p_v_slp = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_slp) ViewFloat3D("pointer_view_slp", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_slp = create_mirror_view(*p_v_slp);
    read_jra(irec, 365, "slp.2016_daily.nc", h_v_slp.data());
    Kokkos::deep_copy(*p_v_slp, h_v_slp);
  } else if (CppParamMod::mytid == 5) {
      p_v_q10 = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_q10) ViewFloat3D("pointer_view_q10", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_q10 = create_mirror_view(*p_v_q10);
    read_jra(irec, 365, "q_10.2016_daily.nc", h_v_q10.data());
    Kokkos::deep_copy(*p_v_q10, h_v_q10);
  } else if (CppParamMod::mytid == 6) {
      p_v_swhf = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_swhf) ViewFloat3D("pointer_view_swhf", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_swhf = create_mirror_view(*p_v_swhf);
    read_jra(irec, 365, "rsds.2016_daily.nc", h_v_swhf.data());
    Kokkos::deep_copy(*p_v_swhf, h_v_swhf);
  } else if (CppParamMod::mytid == 7) {
      p_v_lwhf = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_lwhf) ViewFloat3D("pointer_view_lwhf", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_lwhf = create_mirror_view(*p_v_lwhf);
    read_jra(irec, 365, "rlds.2016_daily.nc", h_v_lwhf.data());
    Kokkos::deep_copy(*p_v_lwhf, h_v_lwhf);
  } else if (CppParamMod::mytid == 8) {
      p_v_precr = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_precr) ViewFloat3D("pointer_view_precr", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_precr = create_mirror_view(*p_v_precr);
    read_jra(irec, 365, "rain.2016_daily.nc", h_v_precr.data());
    Kokkos::deep_copy(*p_v_precr, h_v_precr);
  } else if (CppParamMod::mytid == 9) {
      p_v_precs = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_precs) ViewFloat3D("pointer_view_precs", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_precs = create_mirror_view(*p_v_precs);
    read_jra(irec, 365, "snow.2016_daily.nc", h_v_precs.data());
    Kokkos::deep_copy(*p_v_precs, h_v_precs);
  } else if (CppParamMod::mytid == 10) {
      p_v_si = (ViewFloat3D *) malloc(sizeof(ViewFloat3D));
      new (p_v_si) ViewFloat3D("pointer_view_si", 365, S_JMT, S_IMT);
      ViewFloat3D::HostMirror h_v_si = create_mirror_view(*p_v_si);
    read_jra(irec, 365, "ice.2016_daily.nc", h_v_si.data());
    Kokkos::deep_copy(*p_v_si, h_v_si);
  }
#else
  if (CppParamMod::mytid == 0) {
    p_v_rf = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::rf = new float[365 * S_IMT * S_JMT];
    new (p_v_rf) ViewFloat3D(CppForcMod::rf, 365, S_JMT, S_IMT);
    // update s_lon s_lat
    read_jra(irec, 365, "runoff_all.2016_daily.nc", CppForcMod::buffer, 
        (*p_v_s_lon).data(), (*p_v_s_lat).data(), (*p_v_rf).data());
  } else if (CppParamMod::mytid == 1) {
    p_v_t10 = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::t10 = new float[365 * S_IMT * S_JMT];
    new (p_v_t10) ViewFloat3D(CppForcMod::t10, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "t_10.2016_daily.nc", (*p_v_t10).data());
  } else if (CppParamMod::mytid == 2) {
    p_v_u10 = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::u10 = new float[365 * S_IMT * S_JMT];
    new (p_v_u10) ViewFloat3D(CppForcMod::u10, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "u_10.2016_daily.nc", (*p_v_u10).data());
  } else if (CppParamMod::mytid == 3) {
    p_v_v10 = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::v10 = new float[365 * S_IMT * S_JMT];
    new (p_v_v10) ViewFloat3D(CppForcMod::v10, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "v_10.2016_daily.nc", (*p_v_v10).data());
  } else if (CppParamMod::mytid == 4) {
    p_v_slp = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::slp = new float[365 * S_IMT * S_JMT];
    new (p_v_slp) ViewFloat3D(CppForcMod::slp, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "slp.2016_daily.nc", (*p_v_slp).data());
  } else if (CppParamMod::mytid == 5) {
    p_v_q10 = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::q10 = new float[365 * S_IMT * S_JMT];
    new (p_v_q10) ViewFloat3D(CppForcMod::q10, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "q_10.2016_daily.nc", (*p_v_q10).data());
  } else if (CppParamMod::mytid == 6) {
    p_v_swhf = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::swhf = new float[365 * S_IMT * S_JMT];
    new (p_v_swhf) ViewFloat3D(CppForcMod::swhf, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "rsds.2016_daily.nc", (*p_v_swhf).data());
  } else if (CppParamMod::mytid == 7) {
    p_v_lwhf = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::lwhf = new float[365 * S_IMT * S_JMT];
    new (p_v_lwhf) ViewFloat3D(CppForcMod::lwhf, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "rlds.2016_daily.nc", (*p_v_lwhf).data());
  } else if (CppParamMod::mytid == 8) {
    p_v_precr = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::precr = new float[365 * S_IMT * S_JMT];
    new (p_v_precr) ViewFloat3D(CppForcMod::precr, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "rain.2016_daily.nc", (*p_v_precr).data());
  } else if (CppParamMod::mytid == 9) {
    p_v_precs = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::precs = new float[365 * S_IMT * S_JMT];
    new (p_v_precs) ViewFloat3D(CppForcMod::precs, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "snow.2016_daily.nc", (*p_v_precs).data());
  } else if (CppParamMod::mytid == 10) {
    p_v_si = (ViewFloat3D *) malloc (sizeof (ViewFloat3D));
    CppForcMod::si = new float[365 * S_IMT * S_JMT];
    new (p_v_si) ViewFloat3D(CppForcMod::si, 365, S_JMT, S_IMT);
    read_jra(irec, 365, "ice.2016_daily.nc", (*p_v_si).data());
  }
#endif
  // Broadcasting the s_lon and s_lat
  MPI_Request request_lon, request_lat;

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  MPI_Ibcast(h_v_s_lon.data(), S_IMT * S_JMT, MPI_DOUBLE, 0,
      CppDomain::POP_haloClinic_C.communicator, &request_lon);
  MPI_Ibcast(h_v_s_lat.data(), S_IMT * S_JMT, MPI_DOUBLE, 0,
    CppDomain::POP_haloClinic_C.communicator, &request_lat);
#else
  MPI_Ibcast((*p_v_s_lon).data(), S_IMT * S_JMT, MPI_DOUBLE, 0,
      CppDomain::POP_haloClinic_C.communicator, &request_lon);
  MPI_Ibcast((*p_v_s_lat).data(), S_IMT * S_JMT, MPI_DOUBLE, 0,
    CppDomain::POP_haloClinic_C.communicator, &request_lat);
#endif

  MPI_Wait(&request_lon, MPI_STATUS_IGNORE);
  MPI_Wait(&request_lat, MPI_STATUS_IGNORE);

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  Kokkos::deep_copy(*p_v_s_lon, h_v_s_lon);
  Kokkos::deep_copy(*p_v_s_lat, h_v_s_lat);
#endif
  return ;
}

void kokkos_jra_daily_low (const int &iday) {

  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;

  using KokkosPconstMod::p_v_s_lon;
  using KokkosPconstMod::p_v_s_lat;

  using CppForcMod::t10;
  using CppForcMod::u10;
  using CppForcMod::v10;
  using CppForcMod::slp;
  using CppForcMod::q10;
  using CppForcMod::swhf;
  using CppForcMod::lwhf;
  using CppForcMod::precr;
  using CppForcMod::precs;
  using CppForcMod::rf;
  using CppForcMod::si;

  using CppForcMod::buffer;

  using CppForcMod::tsa3;
  using CppForcMod::wspdu3;
  using CppForcMod::wspdv3;
  using CppForcMod::psa3;
  using CppForcMod::qar3;
  using CppForcMod::swv3;
  using CppForcMod::lwv3;
  using CppForcMod::rain3;
  using CppForcMod::snow3;
  using CppForcMod::runoff3;
  using CppForcMod::seaice3;

  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::JST;
  using CppParamMod::S_IMT;
  using CppParamMod::S_JMT;

  using KokkosForcMod::p_v_t10;
  using KokkosForcMod::p_v_u10;
  using KokkosForcMod::p_v_v10;
  using KokkosForcMod::p_v_slp;
  using KokkosForcMod::p_v_q10;
  using KokkosForcMod::p_v_swhf;
  using KokkosForcMod::p_v_lwhf;
  using KokkosForcMod::p_v_precr;
  using KokkosForcMod::p_v_precs;
  using KokkosForcMod::p_v_rf;
  using KokkosForcMod::p_v_si;
  using KokkosForcMod::p_v_tsa3;
  using KokkosForcMod::p_v_psa3;
  using KokkosForcMod::p_v_qar3;
  using KokkosForcMod::p_v_swv3;
  using KokkosForcMod::p_v_lwv3;
  using KokkosForcMod::p_v_rain3;
  using KokkosForcMod::p_v_snow3;
  using KokkosForcMod::p_v_runoff3;
  using KokkosForcMod::p_v_seaice3;
  using KokkosForcMod::p_v_wspd3;
  using KokkosForcMod::p_v_wspdu3;
  using KokkosForcMod::p_v_wspdv3;
  using KokkosForcMod::p_v_zz;
  using KokkosForcMod::p_v_qs;
  using KokkosForcMod::p_v_ustar;
  using KokkosForcMod::p_v_theta;
  using KokkosForcMod::p_v_core_tau;
  using KokkosForcMod::p_v_core_latent;
  using KokkosForcMod::p_v_core_sensible;
  using KokkosForcMod::p_v_model_sst;
  using KokkosForcMod::p_v_buffer;
  using KokkosPconstMod::p_v_vit;

  int irec = iday;
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  static ViewDouble1D::HostMirror h_v_buffer = create_mirror_view(*p_v_buffer);

  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_tsa3 (tsa3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_wspdu3 (wspdu3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_wspdv3 (wspdv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_psa3 (psa3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_qar3 (qar3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_swv3 (swv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_lwv3 (lwv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_rain3 (rain3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_snow3 (snow3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_runoff3 (runoff3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_seaice3 (seaice3, JMT, IMT);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    printf("jra daily, irec = %d\n", irec);

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_t10 (t10, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_u10 (u10, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_v10 (v10, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_slp (slp, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_q10 (q10, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_swhf (swhf, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_lwhf (lwhf, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_precr (precr, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_precs (precs, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_si (si, 1, S_JMT, S_IMT);
    static Kokkos::View<float ***, Layout, 
        Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
            h_v_rf (rf, 1, S_JMT, S_IMT);

    // ! read in jra data
    read_jra(irec, 1, "t_10.2016_daily.nc", h_v_t10.data());
    Kokkos::deep_copy(*p_v_t10, h_v_t10);
    read_jra(irec, 1, "u_10.2016_daily.nc", h_v_u10.data());
    Kokkos::deep_copy(*p_v_u10, h_v_u10);
    read_jra(irec, 1, "v_10.2016_daily.nc", h_v_v10.data());
    Kokkos::deep_copy(*p_v_v10, h_v_v10);
    read_jra(irec, 1, "slp.2016_daily.nc",  h_v_slp.data());
    Kokkos::deep_copy(*p_v_slp, h_v_slp);
    read_jra(irec, 1, "q_10.2016_daily.nc", h_v_q10.data());
    Kokkos::deep_copy(*p_v_q10, h_v_q10);
    read_jra(irec, 1, "rsds.2016_daily.nc", h_v_swhf.data());
    Kokkos::deep_copy(*p_v_swhf, h_v_swhf);
    read_jra(irec, 1, "rlds.2016_daily.nc", h_v_lwhf.data());
    Kokkos::deep_copy(*p_v_lwhf, h_v_lwhf);
    read_jra(irec, 1, "rain.2016_daily.nc", h_v_precr.data());
    Kokkos::deep_copy(*p_v_precr, h_v_precr);
    read_jra(irec, 1, "snow.2016_daily.nc", h_v_precs.data());
    Kokkos::deep_copy(*p_v_precs, h_v_precs);
    read_jra(irec, 1, "ice.2016_daily.nc",  h_v_si.data());
    Kokkos::deep_copy(*p_v_si, h_v_si);
    read_jra(irec, 1, "runoff_all.2016_daily.nc", h_v_rf.data());
    Kokkos::deep_copy(*p_v_rf, h_v_rf);
#else
    read_jra(irec, 1, "t_10.2016_daily.nc", (*p_v_t10).data());
    read_jra(irec, 1, "u_10.2016_daily.nc", (*p_v_u10).data());
    read_jra(irec, 1, "v_10.2016_daily.nc", (*p_v_v10).data());
    read_jra(irec, 1, "slp.2016_daily.nc",  (*p_v_slp).data());
    read_jra(irec, 1, "q_10.2016_daily.nc", (*p_v_q10).data());
    read_jra(irec, 1, "rsds.2016_daily.nc", (*p_v_swhf).data());
    read_jra(irec, 1, "rlds.2016_daily.nc", (*p_v_lwhf).data());
    read_jra(irec, 1, "rain.2016_daily.nc", (*p_v_precr).data());
    read_jra(irec, 1, "snow.2016_daily.nc", (*p_v_precs).data());
    read_jra(irec, 1, "ice.2016_daily.nc",  (*p_v_si).data());
    read_jra(irec, 1, "runoff_all.2016_daily.nc", (*p_v_rf).data());
#endif
  }

  // interplate to T grid
  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_t10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_t10)));

    parallel_for ("near_t10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_t10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_tsa3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_tsa3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_tsa3, h_v_tsa3);
#else
  scatter_global_jra_r8_((*p_v_tsa3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_tsa3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_u10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_u10)));
    parallel_for ("near_u10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_u10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_wspdu3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_wspdu3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_wspdu3, h_v_wspdu3);
#else
  scatter_global_jra_r8_((*p_v_wspdu3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_wspdu3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_v10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_v10)));
    parallel_for ("near_v10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_v10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_wspdv3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_wspdv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_wspdv3, h_v_wspdv3);
#else
  scatter_global_jra_r8_((*p_v_wspdv3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_wspdv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_slp", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_slp)));
    parallel_for ("near_slp_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_slp_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_psa3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_psa3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_psa3, h_v_psa3);
#else
  scatter_global_jra_r8_((*p_v_psa3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_psa3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_q10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_q10)));
    parallel_for ("near_q10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_q10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_qar3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_qar3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_qar3, h_v_qar3);
#else
  scatter_global_jra_r8_((*p_v_qar3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_qar3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_swhf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_swhf)));
    parallel_for ("near_swhf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_swhf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_swv3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_swv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_swv3, h_v_swv3);
#else
  scatter_global_jra_r8_((*p_v_swv3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_swv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_lwhf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_lwhf)));
    parallel_for ("near_lwhf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_lwhf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_lwv3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_lwv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_lwv3, h_v_lwv3);
#else
  scatter_global_jra_r8_((*p_v_lwv3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_lwv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_precr", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_precr)));
    parallel_for ("near_precr_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_precr_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_rain3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_rain3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_rain3, h_v_rain3);
#else
  scatter_global_jra_r8_((*p_v_rain3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_rain3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_precs", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_precs)));
    parallel_for ("near_precs_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_precs_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_snow3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_snow3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_snow3, h_v_snow3);
#else
  scatter_global_jra_r8_((*p_v_snow3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_snow3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_rf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_rf)));
    parallel_for ("near_rf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_rf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_runoff3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_runoff3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_runoff3, h_v_runoff3);
#else
  scatter_global_jra_r8_((*p_v_runoff3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_runoff3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    parallel_for ("interplation_nearest_si", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(1, (*p_v_si)));
    parallel_for ("near_si_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_si_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_seaice3.data(), h_v_buffer.data(), &CppParamMod::master_task);
  pop_halo_update(h_v_seaice3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_seaice3, h_v_seaice3);
#else
  scatter_global_jra_r8_((*p_v_seaice3).data(), (*p_v_buffer).data(), &CppParamMod::master_task);
  pop_halo_update((*p_v_seaice3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif
  
  parallel_for ("jra_daily_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorJRADaily1());

  parallel_for ("ncar_ocean_fluxes_jra", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorNcarOceanFluxesJra(
          *p_v_wspd3, *p_v_theta, *p_v_model_sst, *p_v_qar3, *p_v_qs, *p_v_zz,
              *p_v_vit, *p_v_core_sensible, *p_v_core_latent, *p_v_core_tau,
                  *p_v_ustar));

  parallel_for ("jra_daily_2", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorJRADaily2());


  // ! tau to U/V grid
  parallel_for ("jra_daily_3", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{JST-1, 1}, koArr2D{JEM, IMM}, tile2D), FunctorJRADaily3());
      // koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorJRADaily3());
  return ;
}

void kokkos_jra_daily_high(const int &iday) {
  using Kokkos::parallel_for;
  using Kokkos::MDRangePolicy;
  using koArr2D = Kokkos::Array<int64_t, 2>;
  using koArr3D = Kokkos::Array<int64_t, 3>;

  using CppForcMod::buffer;

  using CppForcMod::tsa3;
  using CppForcMod::wspdu3;
  using CppForcMod::wspdv3;
  using CppForcMod::psa3;
  using CppForcMod::qar3;
  using CppForcMod::swv3;
  using CppForcMod::lwv3;
  using CppForcMod::rain3;
  using CppForcMod::snow3;
  using CppForcMod::runoff3;
  using CppForcMod::seaice3;

  using CppParamMod::IMT;
  using CppParamMod::JMT;
  using CppParamMod::JST;
  using CppParamMod::S_IMT;
  using CppParamMod::S_JMT;

  using KokkosForcMod::p_v_t10;
  using KokkosForcMod::p_v_u10;
  using KokkosForcMod::p_v_v10;
  using KokkosForcMod::p_v_slp;
  using KokkosForcMod::p_v_q10;
  using KokkosForcMod::p_v_swhf;
  using KokkosForcMod::p_v_lwhf;
  using KokkosForcMod::p_v_precr;
  using KokkosForcMod::p_v_precs;
  using KokkosForcMod::p_v_rf;
  using KokkosForcMod::p_v_si;
  using KokkosForcMod::p_v_tsa3;
  using KokkosForcMod::p_v_psa3;
  using KokkosForcMod::p_v_qar3;
  using KokkosForcMod::p_v_swv3;
  using KokkosForcMod::p_v_lwv3;
  using KokkosForcMod::p_v_rain3;
  using KokkosForcMod::p_v_snow3;
  using KokkosForcMod::p_v_runoff3;
  using KokkosForcMod::p_v_seaice3;
  using KokkosForcMod::p_v_wspd3;
  using KokkosForcMod::p_v_wspdu3;
  using KokkosForcMod::p_v_wspdv3;

  using KokkosForcMod::p_v_zz;
  using KokkosForcMod::p_v_qs;
  using KokkosForcMod::p_v_ustar;
  using KokkosForcMod::p_v_theta;
  using KokkosForcMod::p_v_core_tau;
  using KokkosForcMod::p_v_model_sst;
  using KokkosForcMod::p_v_core_latent;
  using KokkosForcMod::p_v_core_sensible;
  using KokkosForcMod::p_v_buffer;
  using KokkosPconstMod::p_v_vit;
  using CppPOPHaloMod::pop_halo_update_2dr8;

  int irec = iday;
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  static ViewDouble1D::HostMirror h_v_buffer = create_mirror_view(*p_v_buffer);

  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_tsa3 (tsa3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_wspdu3 (wspdu3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_wspdv3 (wspdv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_psa3 (psa3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_qar3 (qar3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_swv3 (swv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_lwv3 (lwv3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_rain3 (rain3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_snow3 (snow3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_runoff3 (runoff3, JMT, IMT);
  static Kokkos::View<double **, Layout, 
      Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> 
          h_v_seaice3 (seaice3, JMT, IMT);
#endif

  if (CppParamMod::mytid == CppParamMod::master_task) {
    printf("jra daily, irec = %d\n", irec);
  }

  int scatter_tid;
  // interplate to T grid
  // rf
  if (CppParamMod::mytid == 0) {
    parallel_for ("interplation_nearest_rf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_rf)));
    parallel_for ("near_rf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_rf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 0; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_runoff3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_runoff3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_runoff3, h_v_runoff3);
#else
  scatter_global_jra_r8_((*p_v_runoff3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_runoff3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // t10
  if (CppParamMod::mytid == 1) {
    parallel_for ("interplation_nearest_t10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_t10)));
    parallel_for ("near_t10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_t10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 1; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_tsa3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_tsa3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_tsa3, h_v_tsa3);
#else
  scatter_global_jra_r8_((*p_v_tsa3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_tsa3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // u10
  if (CppParamMod::mytid == 2) {
    parallel_for ("interplation_nearest_u10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_u10)));
    parallel_for ("near_u10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_u10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 2; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_wspdu3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_wspdu3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_wspdu3, h_v_wspdu3);
#else
  scatter_global_jra_r8_((*p_v_wspdu3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_wspdu3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // v10
  if (CppParamMod::mytid == 3) {
    parallel_for ("interplation_nearest_v10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_v10)));
    parallel_for ("near_v10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_v10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 3; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_wspdv3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_wspdv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_wspdv3, h_v_wspdv3);
#else
  scatter_global_jra_r8_((*p_v_wspdv3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_wspdv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // slp
  if (CppParamMod::mytid == 4) {
    parallel_for ("interplation_nearest_slp", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_slp)));
    parallel_for ("near_slp_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_slp_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 4; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_psa3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_psa3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_psa3, h_v_psa3);
#else
  scatter_global_jra_r8_((*p_v_psa3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_psa3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // q10
  if (CppParamMod::mytid == 5) {
    parallel_for ("interplation_nearest_q10", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_q10)));
    parallel_for ("near_q10_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_q10_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 5; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_qar3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_qar3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_qar3, h_v_qar3);
#else
  scatter_global_jra_r8_((*p_v_qar3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_qar3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // swhf
  if (CppParamMod::mytid == 6) {
    parallel_for ("interplation_nearest_swhf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_swhf)));
    parallel_for ("near_swhf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_swhf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 6; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_swv3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_swv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_swv3, h_v_swv3);
#else
  scatter_global_jra_r8_((*p_v_swv3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_swv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // lwhf
  if (CppParamMod::mytid == 7) {
    parallel_for ("interplation_nearest_lwhf", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_lwhf)));
    parallel_for ("near_lwhf_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_lwhf_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 7; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_lwv3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_lwv3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_lwv3, h_v_lwv3);
#else
  scatter_global_jra_r8_((*p_v_lwv3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_lwv3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // precr
  if (CppParamMod::mytid == 8) {
    parallel_for ("interplation_nearest_precr", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_precr)));
    parallel_for ("near_precr_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_precr_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 8; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_rain3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_rain3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_rain3, h_v_rain3);
#else
  scatter_global_jra_r8_((*p_v_rain3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_rain3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // precs
  if (CppParamMod::mytid == 9) {
    parallel_for ("interplation_nearest_precs", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_precs)));
    parallel_for ("near_precs_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_precs_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 9; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_snow3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_snow3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_snow3, h_v_snow3);
#else
  scatter_global_jra_r8_((*p_v_snow3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_snow3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif

  // si
  if (CppParamMod::mytid == 10) {
    parallel_for ("interplation_nearest_si", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{S_JMT, S_IMT}, tile2D), FunctorInterplationNearest(irec, (*p_v_si)));
    parallel_for ("near_si_1", MDRangePolicy<Kokkos::Rank<2>> (
        koArr2D{0, 0}, koArr2D{JMT_GLOBAL, IMT_GLOBAL}, tile2D), FunctorNear1());
    parallel_for ("near_si_2", JMT_GLOBAL, FunctorNear2());
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
    Kokkos::deep_copy(h_v_buffer, *p_v_buffer);
#endif
  }
  scatter_tid = 10; 
#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
  scatter_global_jra_r8_(h_v_seaice3.data(), h_v_buffer.data(), &scatter_tid);
  pop_halo_update(h_v_seaice3.data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
  Kokkos::deep_copy(*p_v_seaice3, h_v_seaice3);
#else
  scatter_global_jra_r8_((*p_v_seaice3).data(), (*p_v_buffer).data(), &scatter_tid);
  pop_halo_update((*p_v_seaice3).data(), CppParamMod::IMT, CppDomain::POP_haloClinic_C,
      CppPOPGridHorzMod::FLAG_POP_GRID_HORZ_LOC_CENTER,
      CppPOPGridHorzMod::FLAG_POP_FIELD_KIND_SCALAR);
#endif
  
  parallel_for ("jra_daily_1", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorJRADaily1());

  parallel_for ("ncar_ocean_fluxes_jra", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorNcarOceanFluxesJra(
          *p_v_wspd3, *p_v_theta, *p_v_model_sst, *p_v_qar3, *p_v_qs, *p_v_zz,
              *p_v_vit, *p_v_core_sensible, *p_v_core_latent, *p_v_core_tau,
                  *p_v_ustar));

  parallel_for ("jra_daily_2", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{0, 0}, koArr2D{JMT, IMT}, tile2D), FunctorJRADaily2());

  // ! tau to U/V grid
  parallel_for ("jra_daily_3", MDRangePolicy<Kokkos::Rank<2>> (
      koArr2D{JST-1, 1}, koArr2D{JEM, IMM}, tile2D), FunctorJRADaily3());
  return ;
}

#endif // LICOM_ENABLE_KOKKOS

#ifndef LICOM3_KOKKKOS_SRC_UTIL_POP_HALOUPDATE_HPP_
#define LICOM3_KOKKKOS_SRC_UTIL_POP_HALOUPDATE_HPP_

#include "../head/def-undef.h"

#ifdef LICOM_ENABLE_KOKKOS
#include "../head/cpp_param_mod.h"
#include "../head/cpp_dyn_mod.h"
#include "../head/cpp_work_mod.h"
#include "../head/cpp_tracer_mod.h"
#include "../head/cpp_extern_functions.h"

#include "../head/kokkos_config.hpp"

using CppParamMod::IMT;
using CppParamMod::JMT;
using CppParamMod::KM;
using CppParamMod::NTRA;
using CppParamMod::MAX_BLOCKS_CLINIC;

class functor_haloupdate_d2h_i {
 public:
  functor_haloupdate_d2h_i(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble3D &v_array) 
          : arr_stride_k_(0), n_layers_(n_layers), v_buffer_(v_buffer), v_array_3D_(v_array) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int &ii) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_buffer_(ii) = v_array_3D_(iblock, n_layers_ + j, i);
    } else {
      v_buffer_(ii) = v_array_3D_(iblock, JMT - (n_layers_ << 1) - n_layers_ + j, i);
    }
    return ;
  }

  functor_haloupdate_d2h_i(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(0), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array) {}

  functor_haloupdate_d2h_i(const int &arr_stride_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(arr_stride_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int &ii, const int &k) 
      const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = 
          v_array_4D_(iblock, arr_stride_k_ + k, n_layers_ + j, i);
    } else {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = 
          v_array_4D_(iblock, arr_stride_k_ + k, JMT - (n_layers_ << 1) - n_layers_ + j, i);
    }
    return ;
  }
 private:
  const int arr_stride_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble3D v_array_3D_;
  const ViewDouble4D v_array_4D_;
};

class functor_haloupdate_h2d_i {
 public:
  functor_haloupdate_h2d_i(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble3D &v_array) 
          : arr_stride_k_(0), n_layers_(n_layers), v_buffer_(v_buffer), v_array_3D_(v_array) {}
  KOKKOS_INLINE_FUNCTION void operator() (const int &ii) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_array_3D_(iblock, j, i) = v_buffer_(ii);
    } else {
      v_array_3D_(iblock, JMT - (n_layers_ << 1) + j, i) = v_buffer_(ii);
    }
    return ;
  }

  functor_haloupdate_h2d_i(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(0), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array) {}
  functor_haloupdate_h2d_i(const int &arr_stride_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(arr_stride_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array) {}
  KOKKOS_INLINE_FUNCTION void operator() (const int &ii, const int &k) 
      const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_array_4D_(iblock, arr_stride_k_ + k, j, i) = 
          v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    } else {
      v_array_4D_(iblock, arr_stride_k_ + k, JMT - (n_layers_ << 1) + j, i) = 
          v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    }
    return ;
  }

 private:
  const int arr_stride_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble3D v_array_3D_;
  const ViewDouble4D v_array_4D_;
};
class functor_haloupdate_d2h_j {
 public:
  functor_haloupdate_d2h_j(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble3D &v_array) 
          : arr_stride_k_(0), len_k_(1), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_3D_(v_array) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int &jj) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int stride = (n_layers_ * IMT) << 1;
    if (i < n_layers_) {
      v_buffer_(stride + jj) = v_array_3D_(iblock, j, n_layers_ + i);
    } else {
      v_buffer_(stride + jj) = v_array_3D_(iblock, j, IMT - (n_layers_ << 1) - n_layers_ + i);
    }
    return ;
  }

  functor_haloupdate_d2h_j(const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(0), len_k_(len_k), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_4D_(v_array) {}
  functor_haloupdate_d2h_j(const int &arr_stride_k, const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(arr_stride_k), len_k_(len_k), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_4D_(v_array) {}

  KOKKOS_INLINE_FUNCTION void operator() (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int total_layers = n_layers_ << 1;
    const int stride = total_layers * IMT * len_k_;
    if (i < n_layers_) {
      v_buffer_(stride + k * total_layers * JMT + jj) = 
          v_array_4D_(iblock, arr_stride_k_ + k, j, n_layers_ + i);
    } else {
      v_buffer_(stride + k * total_layers * JMT + jj) = 
          v_array_4D_(iblock, arr_stride_k_ + k, j, IMT - total_layers - n_layers_ + i);
    }
    return ;
  }
 private:
  const int arr_stride_k_, len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble3D v_array_3D_;
  const ViewDouble4D v_array_4D_;
};
class functor_haloupdate_h2d_j {
 public:
  functor_haloupdate_h2d_j(const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble3D &v_array) 
          : arr_stride_k_(0), len_k_(1), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_3D_(v_array) {}
  KOKKOS_INLINE_FUNCTION void operator() (const int &jj) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int stride = (n_layers_ * IMT) << 1;
    if (i < n_layers_) {
      v_array_3D_(iblock, j, i) = v_buffer_(stride + jj);
    } else {
      v_array_3D_(iblock, j, IMT - (n_layers_ << 1) + i) = v_buffer_(stride + jj);
    }
    return ;
  }

  functor_haloupdate_h2d_j(const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(0), len_k_(len_k), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_4D_(v_array) {}
  functor_haloupdate_h2d_j(const int &arr_stride_k, const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array) 
          : arr_stride_k_(arr_stride_k), len_k_(len_k), n_layers_(n_layers), 
              v_buffer_(v_buffer), v_array_4D_(v_array) {}
  KOKKOS_INLINE_FUNCTION void operator() (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int total_layers = n_layers_ << 1;
    const int stride = total_layers * IMT * len_k_;
    if (i < n_layers_) {
      v_array_4D_(iblock, arr_stride_k_ + k, j, i) =
          v_buffer_(stride + k * total_layers * JMT + jj);
    } else {
      v_array_4D_(iblock, arr_stride_k_ + k, j, IMT - total_layers + i) =
          v_buffer_(stride + k * total_layers * JMT + jj);
    }
    return ;
  }

 private:
  const int arr_stride_k_, len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble3D v_array_3D_;
  const ViewDouble4D v_array_4D_;
};

class functor_haloupdate_d2h_full_n_layers_i {
 public:
  functor_haloupdate_d2h_full_n_layers_i (const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &ii, const int &k) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = v_array_4D_(iblock, k, j, i);
    } else {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = v_array_4D_(iblock, k, JMT - (n_layers_ << 1) + j, i);
    }
    return ;
  }
 private:
  const int n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};

class functor_haloupdate_d2h_full_n_layers_j {
 public:
  functor_haloupdate_d2h_full_n_layers_j (const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : len_k_(len_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    int total_layers = n_layers_ << 1;
    const int stride = total_layers * len_k_ * IMT;
    if (i < n_layers_) {
      v_buffer_(stride + k * total_layers * JMT + jj) = v_array_4D_(iblock, k, j, i);
    } else {
      v_buffer_(stride + k * total_layers * JMT + jj) = v_array_4D_(iblock, k, j, IMT - total_layers + i);
    }
    return ;
  }
 private:
  const int len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};
class functor_haloupdate_h2d_full_n_layers_i {
 public:
  functor_haloupdate_h2d_full_n_layers_i (const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &ii, const int &k) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_array_4D_(iblock, k, j, i) = v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    } else {
      v_array_4D_(iblock, k, JMT - (n_layers_ << 1) + j, i) = v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    }
    return ;
  }
 private:
  const int n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};
class functor_haloupdate_h2d_full_n_layers_j {
 public:
  functor_haloupdate_h2d_full_n_layers_j (const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : len_k_(len_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int total_layers = n_layers_ << 1;
    const int stride = total_layers * len_k_ * IMT;
    if (i < n_layers_) {
      v_array_4D_(iblock, k, j, i) = v_buffer_(stride + k * total_layers * JMT + jj);
    } else {
      v_array_4D_(iblock, k, j, IMT - total_layers + i) = v_buffer_(stride + k * total_layers * JMT + jj);
    }
    return ;
  }
 private:
  const int len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};

class FuncHaloUpdateD2HTracVtlI {
 public:
  FuncHaloUpdateD2HTracVtlI (const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &ii, const int &k) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = v_array_4D_(iblock, k, j + 2, i);
    } else {
      v_buffer_(k * (n_layers_ << 1) * IMT + ii) = v_array_4D_(iblock, k, JMT - (n_layers_ << 1) - 2 + j, i);
    }
    return ;
  }
 private:
  const int n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};

class FuncHaloUpdateD2HTracVtlJ {
 public:
  FuncHaloUpdateD2HTracVtlJ (const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : len_k_(len_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    int total_layers = n_layers_ << 1;
    const int stride = total_layers * len_k_ * IMT;
    if (i < n_layers_) {
      v_buffer_(stride + k * total_layers * JMT + jj) = v_array_4D_(iblock, k, j, i + 2);
    } else {
      v_buffer_(stride + k * total_layers * JMT + jj) = v_array_4D_(iblock, k, j, IMT - total_layers - 2 + i);
    }
    return ;
  }
 private:
  const int len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};
class FuncHaloUpdateH2DTracVtlI {
 public:
  FuncHaloUpdateH2DTracVtlI (const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &ii, const int &k) const {
    const int iblock = 0;
    const int i = ii % IMT;
    const int j = ii / IMT;
    if (j < n_layers_) {
      v_array_4D_(iblock, k, j, i) = v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    } else {
      v_array_4D_(iblock, k, JMT - (n_layers_ << 1) + j, i) = v_buffer_(k * (n_layers_ << 1) * IMT + ii);
    }
    return ;
  }
 private:
  const int n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};
class FuncHaloUpdateH2DTracVtlJ {
 public:
  FuncHaloUpdateH2DTracVtlJ (const int &len_k, const int &n_layers,
      const ViewDouble1D &v_buffer, const ViewDouble4D &v_array_4D) 
          : len_k_(len_k), n_layers_(n_layers), v_buffer_(v_buffer), v_array_4D_(v_array_4D) {}
  KOKKOS_INLINE_FUNCTION void operator () (const int &jj, const int &k) const {
    const int iblock = 0;
    const int i = jj / JMT;
    const int j = jj % JMT;
    const int total_layers = n_layers_ << 1;
    const int stride = total_layers * len_k_ * IMT;
    if (i < n_layers_) {
      v_array_4D_(iblock, k, j, i) = v_buffer_(stride + k * total_layers * JMT + jj);
    } else {
      v_array_4D_(iblock, k, j, IMT - total_layers + i) = v_buffer_(stride + k * total_layers * JMT + jj);
    }
    return ;
  }
 private:
  const int len_k_, n_layers_;
  const ViewDouble1D v_buffer_;
  const ViewDouble4D v_array_4D_;
};


#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE

using team_policy = Kokkos::TeamPolicy<>;
using TeamHandle = typename team_policy::member_type;

class FuncGetHaloTransDouble {
 public:
  FuncGetHaloTransDouble (const ViewDouble4D &v_src, 
      const ViewDouble1D &v_buff, const int &startLayer, const int &lenLayer,
          const int &lenA, const int &lenB, const int &lenC, const int &maxNumBlock)
    : v_src_(v_src), v_buff_(v_buff), start_layer_(startLayer), len_layer_(lenLayer),
        len_a_(lenA), len_b_(lenB), len_c_(lenC), max_num_block_(maxNumBlock) {}
  KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {

    // using ScratchView = Kokkos::View<double*, 
    //     Kokkos::DefaultExecutionSpace::scratch_memory_space>;

    // team_member.league_rank();
    // team_member.team_rank();
    // team_member.team_barrier();
    const int tid = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    const int idx_a = tid / max_num_block_;

    if (idx_a >= len_a_) {
      return ;
    }
  //   ScratchView v_tile(team_member.team_scratch(0), 
  //       team_member.team_size());

    const int idx_bc = tid % max_num_block_;

    const int len_blk_c0 = len_c_ - (start_layer_ << 1);
    const int len_blk_c1 = len_layer_;

    const int idx_b0 = idx_bc / len_blk_c0;
    const int idx_b1 = idx_bc / len_blk_c1;

    const int idx_c0 = idx_bc % len_blk_c0;
    const int idx_c1 = idx_bc % len_blk_c1;

    // Block 0
    int sb = start_layer_;
    int eb = start_layer_ + len_layer_;
    int sc = start_layer_;
    int ec = len_c_ - start_layer_;

    int num_block = len_layer_ * len_blk_c0;

    if (((sb + idx_b0) < eb) && ((sc + idx_c0) < ec)) {
      v_buff_(idx_a*num_block + idx_b0*len_blk_c0 + idx_c0) = 
          v_src_(0, idx_a, sb+idx_b0, sc+idx_c0);
    }

    // Block 1
    int start_buf = num_block * len_a_;

    sb = start_layer_ + len_layer_;
    eb = len_b_ - start_layer_ - len_layer_;
    sc = start_layer_;
    ec = start_layer_ + len_layer_;

    num_block = (eb - sb) * len_blk_c1;
  
    if (((sb + idx_b1) < eb) && ((sc + idx_c1) < ec)) {
      v_buff_(start_buf + idx_a*num_block + idx_b1*len_blk_c1 + idx_c1) = 
          v_src_(0, idx_a, sb+idx_b1, sc+idx_c1);
    }
 
    // Block 2
    start_buf += (num_block * len_a_);

    sb = start_layer_ + len_layer_;
    eb = len_b_ - start_layer_ - len_layer_;
    sc = len_c_ - start_layer_ - len_layer_;
    ec = len_c_ - start_layer_;

    num_block = (eb - sb) * len_blk_c1;
  
    if (((sb + idx_b1) < eb) && ((sc + idx_c1) < ec)) {
      v_buff_(start_buf + idx_a*num_block + idx_b1*len_blk_c1 + idx_c1) = 
          v_src_(0, idx_a, sb+idx_b1, sc+idx_c1);
    }
 
    // Block 3
    start_buf += (num_block * len_a_);

    sb = len_b_ - start_layer_ - len_layer_;
    eb = len_b_ - start_layer_;
    sc = start_layer_;
    ec = len_c_ - start_layer_;

    num_block = len_layer_ * len_blk_c0;
  
    if (((sb + idx_b0) < eb) && ((sc + idx_c0) < ec)) {
      v_buff_(start_buf + idx_a*num_block + idx_b0*len_blk_c0 + idx_c0) = 
          v_src_(0, idx_a, sb+idx_b0, sc+idx_c0);
    }
    return ;
  }
 private:
  const ViewDouble4D v_src_;
  const ViewDouble1D v_buff_;
  const int start_layer_, len_layer_;
  const int len_a_, len_b_, len_c_;
  const int max_num_block_;
};

// class FuncGetHaloTransDouble {
//  public:
//   FuncGetHaloTransDouble (const ViewDouble4D &v_src, 
//       const ViewDouble1D &v_buff, const int &startLayer, const int &lenLayer,
//           const int &lenA, const int &lenB, const int &lenC, const int &maxNumBlock)
//     : v_src_(v_src), v_buff_(v_buff), start_layer_(startLayer), len_layer_(lenLayer),
//         len_a_(lenA), len_b_(lenB), len_c_(lenC), max_num_block_(maxNumBlock) {}
//   KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {

//     using ScratchView = Kokkos::View<double*, 
//         Kokkos::DefaultExecutionSpace::scratch_memory_space>;

//     // 1056 = 32 * 33
//     ScratchView v_shmem (team_member.team_scratch(0), 1056);
//     // ScratchView v_shmem (team_member.team_scratch(0), 
//     //     team_member.team_size());

//     // team_member.league_rank();
//     // team_member.team_rank();
//     // team_member.team_barrier();

//     const int tid = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

//     int sb, eb, sc, ec;
//     int len_blk_b, len_blk_c;
//     int num_blk;
//     int idx_a, idx_b, idx_c, idx_blk;
//     int start_buff;

//     // // Block 0
//     sb = start_layer_;
//     eb = start_layer_ + len_layer_;
//     sc = start_layer_;
//     ec = len_c_ - start_layer_;

//     len_blk_b = eb - sb;
//     len_blk_c = ec - sc;

//     // Array -> shmem
//     start_buff = 0;
//     num_blk = len_blk_b * len_blk_c;
//     idx_a   = tid     / num_blk;
//     idx_blk = tid     % num_blk;
//     idx_b   = idx_blk / len_blk_c;
//     idx_c   = idx_blk % len_blk_c;
//     if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//       v_shmem(team_member.team_rank()) = v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     }
//     team_member.team_barrier();
//     // shmem -> buffer
//     // num_blk    = len_blk_c * len_a_;
//     // idx_b      = tid     / num_blk;
//     // idx_blk    = tid     % num_blk;
//     // idx_c      = idx_blk / len_a_;
//     // idx_a      = idx_blk % len_a_;

//     idx_blk    = tid     / len_a_;
//     idx_b      = idx_blk / len_blk_c;
//     idx_c      = idx_blk % len_blk_c;

//     idx_a      = tid     % len_a_;

//     if ((idx_a < len_a_) && (idx_b < len_blk_b) && (idx_c < len_blk_c)) {
//       v_buff_(start_buff + idx_b*len_blk_c*len_a_ + idx_c*len_a_ + idx_a) = 
//           v_shmem(team_member.team_rank());
//     }

//     // // Block 1
//     // start_buff += (len_blk_b * len_blk_c * len_a_);

//     // sb = start_layer_ + len_layer_;
//     // eb = len_b_ - start_layer_ - len_layer_;
//     // sc = start_layer_;
//     // ec = start_layer_ + len_layer_;

//     // len_blk_b = eb - sb;
//     // len_blk_c = ec - sc;

//     // // Array -> shmem
//     // num_blk = len_blk_b * len_blk_c;
//     // idx_a   = tid     / num_blk;
//     // idx_blk = tid     % num_blk;
//     // idx_b   = idx_blk / len_blk_c;
//     // idx_c   = idx_blk % len_blk_c;
//     // if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//     //   v_shmem(team_member.team_rank()) = v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     // }
//     // // shmem -> buffer
//     // // num_blk     = len_blk_c * len_a_;
//     // // idx_b       = tid     / num_blk;
//     // // idx_blk     = tid     % num_blk;
//     // // idx_c       = idx_blk / len_a_;
//     // // idx_a       = idx_blk % len_a_;
//     // // start_buff += (num_blk * len_a_);
//     // idx_blk     = tid     / len_a_;
//     // idx_b       = idx_blk / len_blk_c;
//     // idx_c       = idx_blk % len_blk_c;
//     // idx_a       = tid     % len_a_;
//     // if ((idx_a < len_a_) && (idx_b < len_blk_b) && (idx_c < len_blk_c)) {
//     //   v_buff_(start_buff + idx_b*len_blk_c*len_a_ + idx_c*len_a_ + idx_a) = 
//     //       v_shmem(team_member.team_rank());
//     // }

//     // // Block 2
//     // start_buff += (len_blk_b * len_blk_c * len_a_);
//     // sb = start_layer_ + len_layer_;
//     // eb = len_b_ - start_layer_ - len_layer_;
//     // sc = len_c_ - start_layer_ - len_layer_;
//     // ec = len_c_ - start_layer_;

//     // len_blk_b = eb - sb;
//     // len_blk_c = ec - sc;

//     // // Array -> shmem
//     // num_blk = len_blk_b * len_blk_c;
//     // idx_a   = tid     / num_blk;
//     // idx_blk = tid     % num_blk;
//     // idx_b   = idx_blk / len_blk_c;
//     // idx_c   = idx_blk % len_blk_c;
//     // if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//     //   v_shmem(team_member.team_rank()) = v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     // }
//     // // shmem -> buffer
//     // // num_blk     = len_blk_c * len_a_;
//     // // idx_b       = tid     / num_blk;
//     // // idx_blk     = tid     % num_blk;
//     // // idx_c       = idx_blk / len_a_;
//     // // idx_a       = idx_blk % len_a_;
//     // idx_blk     = tid     / len_a_;
//     // idx_b       = idx_blk / len_blk_c;
//     // idx_b       = idx_blk % len_blk_c;
//     // idx_a       = tid     % len_a_;
//     // if ((idx_a < len_a_) && (idx_b < len_blk_b) && (idx_c < len_blk_c)) {
//     //   v_buff_(start_buff + idx_b*len_blk_c*len_a_ + idx_c*len_a_ + idx_a) = 
//     //       v_shmem(team_member.team_rank());
//     // }

//     // // Block 3
//     // start_buff += (len_blk_b * len_blk_c * len_a_);
//     // sb = len_b_ - start_layer_ - len_layer_;
//     // eb = len_b_ - start_layer_;
//     // sc = start_layer_;
//     // ec = len_c_ - start_layer_;

//     // len_blk_b = eb - sb;
//     // len_blk_c = ec - sc;

//     // // Array -> shmem
//     // num_blk = len_blk_b * len_blk_c;
//     // idx_a   = tid     / num_blk;
//     // idx_blk = tid     % num_blk;
//     // idx_b   = idx_blk / len_blk_c;
//     // idx_c   = idx_blk % len_blk_c;
//     // if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//     //   v_shmem(team_member.team_rank()) = v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     // }
//     // // shmem -> buffer
//     // idx_blk     = tid     / len_a_;
//     // idx_b       = idx_blk / len_blk_c;
//     // idx_c       = idx_blk % len_blk_c;
//     // idx_a       = tid     % len_a_;
//     // if ((idx_a < len_a_) && (idx_b < len_blk_b) && (idx_c < len_blk_c)) {
//     //   v_buff_(start_buff + idx_b*len_blk_c*len_a_ + idx_c*len_a_ + idx_a) = 
//     //       v_shmem(team_member.team_rank());
//     // }
//     // int tid = team_member.league_rank() * team_member.team_size() + team_member.team_rank();

//     // if (idx_a >= len_a_) {
//     //   return ;
//     // }
//   //   ScratchView v_tile(team_member.team_scratch(0), 
//   //       team_member.team_size());

//     idx_a = tid / max_num_block_;
//     int idx_bc = tid % max_num_block_;

//     // Block 0
//     sb = start_layer_;
//     eb = start_layer_ + len_layer_;
//     sc = start_layer_;
//     ec = len_c_ - start_layer_;
//     int len_block_c = ec - sc;
//     idx_b = idx_bc / len_block_c;
//     idx_c = idx_bc % len_block_c;

//     int num_block = (eb - sb) * len_block_c;
//     int stride_buff = len_block_c * len_a_;

//     if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//       v_buff_(idx_b*stride_buff + idx_c*len_a_ + idx_a) = 
//           v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     }

//     // Block 1
//     int start_buf = num_block * len_a_;
//     sb = start_layer_ + len_layer_;
//     eb = len_b_ - start_layer_ - len_layer_;
//     sc = start_layer_;
//     ec = start_layer_ + len_layer_;
//     len_block_c = ec - sc;

//     idx_b = idx_bc / len_block_c;
//     idx_c = idx_bc % len_block_c;

//     num_block = (eb - sb) * len_block_c;
//     stride_buff = len_block_c * len_a_;
  
//     if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//       v_buff_(start_buf + idx_b*stride_buff + idx_c*len_a_ + idx_a) = 
//           v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     }
 
//     // Block 2
//     start_buf += (num_block * len_a_);
//     sb = start_layer_ + len_layer_;
//     eb = len_b_ - start_layer_ - len_layer_;
//     sc = len_c_ - start_layer_ - len_layer_;
//     ec = len_c_ - start_layer_;
//     len_block_c = ec - sc;

//     idx_b = idx_bc / len_block_c;
//     idx_c = idx_bc % len_block_c;

//     num_block = (eb - sb) * len_block_c;
//     stride_buff = len_block_c * len_a_;
  
//     if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//       v_buff_(start_buf + idx_b*stride_buff + idx_c*len_a_ + idx_a) = 
//           v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     }
 
//     // Block 3
//     start_buf += (num_block * len_a_);
//     sb = len_b_ - start_layer_ - len_layer_;
//     eb = len_b_ - start_layer_;
//     sc = start_layer_;
//     ec = len_c_ - start_layer_;
//     len_block_c = ec - sc;

//     idx_b = idx_bc / len_block_c;
//     idx_c = idx_bc % len_block_c;

//     num_block = (eb - sb) * len_block_c;
//     stride_buff = len_block_c * len_a_;
  
//     if ((idx_a < len_a_) && ((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
//       v_buff_(start_buf + idx_b*stride_buff + idx_c*len_a_ + idx_a) = 
//           v_src_(0, idx_a, sb+idx_b, sc+idx_c);
//     }
//     return ;
//   }
//  private:
//   const ViewDouble4D v_src_;
//   const ViewDouble1D v_buff_;
//   const int start_layer_, len_layer_;
//   const int len_a_, len_b_, len_c_;
//   const int max_num_block_;
// };

class FuncPutHaloTransDouble {
 public:
  FuncPutHaloTransDouble (const ViewDouble1D &v_buff, 
      const ViewDouble4D &v_dst,  const int &startLayer, const int &lenLayer,
          const int &lenA, const int &lenB, const int &lenC, const int &maxNumBlock)
    : v_buff_(v_buff), v_dst_(v_dst), start_layer_(startLayer), len_layer_(lenLayer),
        len_a_(lenA), len_b_(lenB), len_c_(lenC), max_num_block_(maxNumBlock) {}

  // KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {

  //   using ScratchView = Kokkos::View<double*, 
  //       Kokkos::DefaultExecutionSpace::scratch_memory_space>;

  //   ScratchView v_tile(team_member.team_scratch(0), 
  //       team_member.team_size());

  //   int tid = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
  //   int idx_a = tid / max_num_block_;

  //   if (idx_a >= len_a_) {
  //     return ;
  //   }

  //   Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, max_num_block_), 
  //       [&](const int &idx_bc) {
  //     int sb = start_layer_;
  //     int eb = start_layer_ + len_layer_;
  //     int sc = start_layer_;
  //     int ec = len_c_ - start_layer_;
  //     int len_block_c = ec - sc;

  //     int idx_b = idx_bc / len_block_c;
  //     int idx_c = idx_bc % len_block_c;

  //     int num_block = (eb - sb) * len_block_c;
   
  //     if (((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
  //       v_dst_(0, idx_a, sb+idx_b, sc+idx_c) =
  //           v_buff_(idx_a*num_block + idx_b*len_block_c + idx_c);
  //     }
   
  //     // Block 1
  //     int start_buf = num_block * len_a_;
  //     sb = start_layer_ + len_layer_;
  //     eb = len_b_ - start_layer_ - len_layer_;
  //     sc = start_layer_;
  //     ec = start_layer_ + len_layer_;
  //     len_block_c = ec - sc;

  //     idx_b = idx_bc / len_block_c;
  //     idx_c = idx_bc % len_block_c;

  //     num_block = (eb - sb) * len_block_c;
   
  //     if (((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
  //       v_dst_(0, idx_a, sb+idx_b, sc+idx_c) =
  //           v_buff_(start_buf + idx_a*num_block + idx_b*len_block_c + idx_c);
  //     }
   
  //     // Block 2
  //     start_buf += (num_block * len_a_);
  //     sb = start_layer_ + len_layer_;
  //     eb = len_b_ - start_layer_ - len_layer_;
  //     sc = len_c_ - start_layer_ - len_layer_;
  //     ec = len_c_ - start_layer_;
  //     len_block_c = ec - sc;

  //     idx_b = idx_bc / len_block_c;
  //     idx_c = idx_bc % len_block_c;

  //     num_block = (eb - sb) * len_block_c;
   
  //     if (((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
  //       v_dst_(0, idx_a, sb+idx_b, sc+idx_c) =
  //           v_buff_(start_buf + idx_a*num_block + idx_b*len_block_c + idx_c);
  //     }
   
  //     // // Block 3
  //     start_buf += (num_block * len_a_);
  //     sb = len_b_ - start_layer_ - len_layer_;
  //     eb = len_b_ - start_layer_;
  //     sc = start_layer_;
  //     ec = len_c_ - start_layer_;
  //     len_block_c = ec - sc;

  //     idx_b = idx_bc / len_block_c;
  //     idx_c = idx_bc % len_block_c;

  //     num_block = (eb - sb) * len_block_c;
   
  //     if (((sb + idx_b) < eb) && ((sc + idx_c) < ec)) {
  //       v_dst_(0, idx_a, sb+idx_b, sc+idx_c) =
  //           v_buff_(start_buf + idx_a*num_block + idx_b*len_block_c + idx_c);
  //     }

  //     return ;
  //   });

  //   return ;
  // }

  KOKKOS_INLINE_FUNCTION void operator () (const TeamHandle& team_member) const {

    const int tid = team_member.league_rank() * team_member.team_size() + team_member.team_rank();
    const int idx_a = tid / max_num_block_;

    if (idx_a >= len_a_) {
      return ;
    }

    const int idx_bc = tid % max_num_block_;

    const int len_blk_c0 = len_c_ - (start_layer_ << 1);
    const int len_blk_c1 = len_layer_;

    const int idx_b0 = idx_bc / len_blk_c0;
    const int idx_b1 = idx_bc / len_blk_c1;

    const int idx_c0 = idx_bc % len_blk_c0;
    const int idx_c1 = idx_bc % len_blk_c1;

    // Block 0
    int sb = start_layer_;
    int eb = start_layer_ + len_layer_;
    int sc = start_layer_;
    int ec = len_c_ - start_layer_;

    int num_block = len_layer_ * len_blk_c0;

    if (((sb + idx_b0) < eb) && ((sc + idx_c0) < ec)) {
      v_dst_(0, idx_a, sb+idx_b0, sc+idx_c0) =
          v_buff_(idx_a*num_block + idx_b0*len_blk_c0 + idx_c0);
    }

    // Block 1
    int start_buf = num_block * len_a_;

    sb = start_layer_ + len_layer_;
    eb = len_b_ - start_layer_ - len_layer_;
    sc = start_layer_;
    ec = start_layer_ + len_layer_;

    num_block = (eb - sb) * len_blk_c1;
  
    if (((sb + idx_b1) < eb) && ((sc + idx_c1) < ec)) {
      v_dst_(0, idx_a, sb+idx_b1, sc+idx_c1) =
          v_buff_(start_buf + idx_a*num_block + idx_b1*len_blk_c1 + idx_c1);
    }
 
    // Block 2
    start_buf += (num_block * len_a_);

    sb = start_layer_ + len_layer_;
    eb = len_b_ - start_layer_ - len_layer_;
    sc = len_c_ - start_layer_ - len_layer_;
    ec = len_c_ - start_layer_;

    num_block = (eb - sb) * len_blk_c1;

    if (((sb + idx_b1) < eb) && ((sc + idx_c1) < ec)) {
      v_dst_(0, idx_a, sb+idx_b1, sc+idx_c1) =
          v_buff_(start_buf + idx_a*num_block + idx_b1*len_blk_c1 + idx_c1);
    }
 
    // Block 3
    start_buf += (num_block * len_a_);

    sb = len_b_ - start_layer_ - len_layer_;
    eb = len_b_ - start_layer_;
    sc = start_layer_;
    ec = len_c_ - start_layer_;

    num_block = len_layer_ * len_blk_c0;
  
    if (((sb + idx_b0) < eb) && ((sc + idx_c0) < ec)) {
      v_dst_(0, idx_a, sb+idx_b0, sc+idx_c0) =
          v_buff_(start_buf + idx_a*num_block + idx_b0*len_blk_c0 + idx_c0);
    }

    return ;
  }

 private:
  const ViewDouble1D v_buff_;
  const ViewDouble4D v_dst_;
  const int start_layer_, len_layer_;
  const int len_a_, len_b_, len_c_;
  const int max_num_block_;
};

#endif // KOKKOS_ENABLE_DEVICE_MEM_SPACE

// KOKKOS_REGISTER_FOR_1D(functor_haloupdate_d2h_i,  functor_haloupdate_d2h_i)
// KOKKOS_REGISTER_FOR_1D(functor_haloupdate_d2h_j,  functor_haloupdate_d2h_j)
// KOKKOS_REGISTER_FOR_1D(functor_haloupdate_h2d_i,  functor_haloupdate_h2d_i)
// KOKKOS_REGISTER_FOR_1D(functor_haloupdate_h2d_j,  functor_haloupdate_h2d_j)

// KOKKOS_REGISTER_FOR_2D(functor_haloupdate_d2h_i,  functor_haloupdate_d2h_i)
// KOKKOS_REGISTER_FOR_2D(functor_haloupdate_d2h_j,  functor_haloupdate_d2h_j)
// KOKKOS_REGISTER_FOR_2D(functor_haloupdate_h2d_i,  functor_haloupdate_h2d_i)
// KOKKOS_REGISTER_FOR_2D(functor_haloupdate_h2d_j,  functor_haloupdate_h2d_j)

// KOKKOS_REGISTER_FOR_2D(FuncHaloUpdateD2HTracVtlI, FuncHaloUpdateD2HTracVtlI)
// KOKKOS_REGISTER_FOR_2D(FuncHaloUpdateD2HTracVtlJ, FuncHaloUpdateD2HTracVtlJ)
// KOKKOS_REGISTER_FOR_2D(FuncHaloUpdateH2DTracVtlI, FuncHaloUpdateH2DTracVtlI)
// KOKKOS_REGISTER_FOR_2D(FuncHaloUpdateH2DTracVtlJ, FuncHaloUpdateH2DTracVtlJ)
#endif // LICOM_ENABLE_KOKKOS

#endif // LICOM3_KOKKKOS_SRC_UTIL_POP_HALOUPDATE_HPP_

#ifndef KOKKOS_ATHREAD_KOKKOS_REGISTERFUNCTION_HPP_
#define KOKKOS_ATHREAD_KOKKOS_REGISTERFUNCTION_HPP_

#include "Kokkos_Athread_ParamWrap.h"

#include "slave.h"

#include <cstdio>
#include <cstdlib>

struct AthradRegisterFunctionListNode {
  int   *key;
  void (*fp)(AthreadParamWrap *);
  AthradRegisterFunctionListNode *next;
  int num_intv16;
};

extern int* g_athread_functor_key __attribute__ ((aligned(64)));
template <typename FunctorType>
void athread_get_key(const int dim, int &num_intv16) {
  const char* full_name = __PRETTY_FUNCTION__;
  // int start = 0;
  // for (int i = 0; full_name[i] != '\n'; ++i) {
  //   if (full_name[i] == '=') {
  //     start = i + 2;
  //     break;
  //   }
  // }
  // printf("full_name: %s, start = %d, %c\n", full_name, start, full_name[start]);
  const int start = 52;
  int key_len = 0;
  for (int i = start + 1; full_name[i] != '\n'; ++i) {
    if (full_name[i] == ']') {
      key_len = i - start;
      break;
    }
  }

  if ((key_len == 0) || (key_len >= 32)) {
    printf("error in athread_get_key_simd: %s, len = %d\n", __PRETTY_FUNCTION__, key_len);
    exit(EXIT_FAILURE);
  }

  for (int i = 0; i < key_len; ++i) {
    g_athread_functor_key[i] = full_name[start + i];
  }
  g_athread_functor_key[key_len] = '0' + dim;
  key_len += 1;
  if (key_len > 16) {
    for (int i = key_len; i < 32; ++i) {
      g_athread_functor_key[i] = '0';
    }
    num_intv16 = 2;
  } else {
    for (int i = key_len; i < 16; ++i) {
      g_athread_functor_key[i] = '0';
    }
    num_intv16 = 1;
  }
  return ;
}

#define KOKKOS_ATHREAD_REGISTER_FOR_1D(FUNC_NAME, ...)                     \
void FUNC_NAME##_1D(AthreadParamWrap* para) {                              \
  const int my_athread_tid      = para->my_athread_tid;                    \
  const int ATHREAD_SLAVE_CORES = para->num_athread_cores;                 \
  const int start = para->range[0][0];                                     \
  const int end   = para->range[0][1];                                     \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);    \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                  \
  const int len   = end - start;                                           \
  const int times = (len + ATHREAD_SLAVE_CORES - 1) / ATHREAD_SLAVE_CORES; \
  for (int i = 0; i < times; ++i) {                                        \
    int index = start + athread_tid * times + i;                           \
    if (index < end) {                                                     \
      F->operator() (index);                                               \
    }                                                                      \
  }                                                                        \
  return ;                                                                 \
}                                                                          \

#define KOKKOS_ATHREAD_REGISTER_FOR_2D(FUNC_NAME, ...)                  \
void FUNC_NAME##_2D(AthreadParamWrap* para) {                           \
  const int my_athread_tid = para->my_athread_tid;                      \
  const int ATHREAD_SLAVE_CORES = para->num_athread_cores;              \
  const int start0 = para->range[0][0];                                 \
  const int end0   = para->range[0][1];                                 \
  const int start1 = para->range[1][0];                                 \
  const int end1   = para->range[1][1];                                 \
  const int tile0  = para->tile[0];                                     \
  const int tile1  = para->tile[1];                                     \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor); \
  auto* F = const_cast<__VA_ARGS__ *>(f);                               \
  const int num_tiles0 = (end0 - start0 + tile0 - 1) / tile0;           \
  const int num_tiles1 = (end1 - start1 + tile1 - 1) / tile1;           \
  const int num_tiles = num_tiles0 * num_tiles1;                        \
  for (int index_tile = my_athread_tid; index_tile < num_tiles;         \
      index_tile += ATHREAD_SLAVE_CORES) {                              \
    const int index_tile0 = index_tile / num_tiles1;                    \
    const int index_tile1 = index_tile % num_tiles1;                    \
    for (int index00 = 0; index00 < tile0; ++ index00) {                \
      const int index0 = start0 + index_tile0 * tile0 + index00;        \
      if (index0 >= end0) { break; }                                    \
      for (int index_11 = 0; index_11 < tile1; ++ index_11) {           \
        const int index1 = start1 + index_tile1 * tile1 + index_11;     \
        if (index1 >= end1) { break; }                                  \
        F->operator() (index0, index1);                                 \
      }                                                                 \
    }                                                                   \
  }                                                                     \
  return;                                                               \
}                                                                       \

#define KOKKOS_ATHREAD_REGISTER_FOR_3D(FUNC_NAME, ...)                            \
void FUNC_NAME##_3D(AthreadParamWrap* para) {                                     \
  const int my_athread_tid      = para->my_athread_tid;                           \
  const int ATHREAD_SLAVE_CORES = para->num_athread_cores;                        \
  const int start[3] = {para->range[0][0], para->range[1][0], para->range[2][0]}; \
  const int end[3]   = {para->range[0][1], para->range[1][1], para->range[2][1]}; \
  const int tile[3]  = {para->tile[0],     para->tile[1],     para->tile[2]};     \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);           \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                         \
  const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];             \
  const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];             \
  const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];             \
  const int tmp0 = num_tiles1 * num_tiles2;                                       \
  const int num_tiles = num_tiles0 * tmp0;                                        \
  for (int index_tile = my_athread_tid; index_tile < num_tiles;                   \
      index_tile += ATHREAD_SLAVE_CORES) {                                        \
    const int index_tile0 = index_tile / tmp0;                                    \
    const int tmp1        = index_tile % tmp0;                                    \
    const int index_tile1 = tmp1      / num_tiles2;                               \
    const int index_tile2 = tmp1      % num_tiles2;                               \
    for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {                     \
      const int index0 = start[0] + index_tile0 * tile[0] + index_00;             \
      if (index0 >= end[0]) { break; }                                            \
      for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {                   \
        const int index1 = start[1] + index_tile1 * tile[1] + index_11;           \
        if (index1 >= end[1]) { break; }                                          \
        for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {                 \
          const int index2 = start[2] + index_tile2 * tile[2] + index_22;         \
          if (index2 >= end[2]) { break; }                                        \
          F->operator() (index0, index1, index2);                                 \
        }                                                                         \
      }                                                                           \
    }                                                                             \
  }                                                                               \
  return ;                                                                        \
}                                                                                 \

#define KOKKOS_ATHREAD_REGISTER_FOR_4D(FUNC_NAME, ...)                      \
void FUNC_NAME##_4D(AthreadParamWrap* para) {                               \
  const int my_athread_tid      = para->my_athread_tid;                     \
  const int ATHREAD_SLAVE_CORES = para->num_athread_cores;                  \
  const int start[4] = {para->range[0][0], para->range[1][0],               \
                        para->range[2][0], para->range[3][0]};              \
  const int end[4]   = {para->range[0][1], para->range[1][1],               \
                        para->range[2][1], para->range[3][1]};              \
  const int tile[4]  = {para->tile[0],     para->tile[1],                   \
                        para->tile[2]    , para->tile[3]};                  \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);     \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                   \
  const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];       \
  const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];       \
  const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];       \
  const int num_tiles3 = (end[3] - start[3] + tile[3] - 1) / tile[3];       \
  const int tmp0       = num_tiles2 * num_tiles3;                           \
  const int tmp1       = num_tiles1 * tmp0;                                 \
  const int num_tiles  = num_tiles0 * tmp1;                                 \
  for (int index_tile = my_athread_tid; index_tile < num_tiles;             \
      index_tile += ATHREAD_SLAVE_CORES) {                                  \
    const int index_tile0 = index_tile / tmp1;                              \
    const int tmp2        = index_tile % tmp1;                              \
    const int index_tile1 = tmp2       / tmp0;                              \
    const int tmp3        = tmp2       % tmp0;                              \
    const int index_tile2 = tmp3       / num_tiles3;                        \
    const int index_tile3 = tmp3       % num_tiles3;                        \
    for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {               \
      const int index0 = start[0] + index_tile0 * tile[0] + index_00;       \
      if (index0 >= end[0]) { break; }                                      \
      for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {             \
        const int index1 = start[1] + index_tile1 * tile[1] + index_11;     \
        if (index1 >= end[1]) { break; }                                    \
        for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {           \
          const int index2 = start[2] + index_tile2 * tile[2] + index_22;   \
          if (index2 >= end[2]) { break; }                                  \
          for (int index_33 = 0; index_33 < tile[3]; ++ index_33) {         \
            const int index3 = start[3] + index_tile3 * tile[3] + index_33; \
            if (index3 >= end[3]) { break; }                                \
            F->operator() (index0, index1, index2, index3);                 \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  return ;                                                                  \
}                                                                           \
                                                                                   
// #define KOKKOS_ATHREAD_REGISTER_FOR_5D(FUNC_NAME, ...)                                  \
// void FUNC_NAME##_5D(AthreadParamWrap* para) {                                           \
//   const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
//   auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
//   F->operator() (                                                                       \
//     para->index[0],                                                                     \
//     para->index[1],                                                                     \
//     para->index[2],                                                                     \
//     para->index[3],                                                                     \
//     para->index[4]                                                                      \
//   );                                                                                    \
// }                                                                                       \

// #define KOKKOS_ATHREAD_REGISTER_FOR_6D(FUNC_NAME, ...)                                  \
// void FUNC_NAME##_6D(AthreadParamWrap* para) {                                           \
//   const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
//   auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
//   F->operator() (                                                                       \
//     para->index[0],                                                                     \
//     para->index[1],                                                                     \
//     para->index[2],                                                                     \
//     para->index[3],                                                                     \
//     para->index[4],                                                                     \
//     para->index[5]                                                                      \
//   );                                                                                    \
// }                                                                                       \

// #define KOKKOS_ATHREAD_REGISTER_REDUCE_1D(FUNC_NAME, ...)                               \
// void FUNC_NAME##_1D(AthreadParamWrap* para) {                                           \
//   const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
//   auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
//   F->operator() (                                                                       \
//     para->index[0],                                                                     \
//     para->reduce_result                                                                 \
//   );                                                                                    \
// }                                                                                       \
                                                                                       
// #define KOKKOS_ATHREAD_REGISTER_REDUCE_2D(FUNC_NAME, ...)                               \
// void FUNC_NAME##_2D(AthreadParamWrap* para) {                                           \
//   const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
//   auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
//   F->operator() (                                                                       \
//     para->index[0],                                                                     \
//     para->index[1],                                                                     \
//     para->reduce_result                                                                 \
//   );                                                                                    \
// }                                                                                       \
                                                                                       
// #define KOKKOS_ATHREAD_REGISTER_REDUCE_3D(FUNC_NAME, ...)                               \
// void FUNC_NAME##_3D(AthreadParamWrap* para) {                                           \
//   const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
//   auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
//   F->operator() (                                                                       \
//     para->index[0],                                                                     \
//     para->index[1],                                                                     \
//     para->index[2],                                                                     \
//     para->reduce_result                                                                 \
//   );                                                                                    \
// }                                                                                       \

#endif // KOKKOS_ATHREAD_REGISTERMACROS_HPP_
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
  // printf("full_name = %s\n", full_name);
  int start = 0;
  for (int i = 0; full_name[i] != '\n'; ++i) {
    if (full_name[i] == '=') {
      start = i + 2;
      break;
    }
  }
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

#define KOKKOS_ATHREAD_REGISTER_FOR_1D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_1D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (para->index[0]);                                                       \
}                                                                                       \

#define KOKKOS_ATHREAD_REGISTER_FOR_2D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_2D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1]                                                                      \
  );                                                                                    \
}                                                                                       \
                                                                               
#define KOKKOS_ATHREAD_REGISTER_FOR_3D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_3D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->index[2]                                                                      \
  );                                                                                    \
}                                                                                       \
                                                                                 
#define KOKKOS_ATHREAD_REGISTER_FOR_4D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_4D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->index[2],                                                                     \
    para->index[3]                                                                      \
  );                                                                                    \
}                                                                                       \
                                                                                   
#define KOKKOS_ATHREAD_REGISTER_FOR_5D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_5D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->index[2],                                                                     \
    para->index[3],                                                                     \
    para->index[4]                                                                      \
  );                                                                                    \
}                                                                                       \

#define KOKKOS_ATHREAD_REGISTER_FOR_6D(FUNC_NAME, ...)                                  \
void FUNC_NAME##_6D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->index[2],                                                                     \
    para->index[3],                                                                     \
    para->index[4],                                                                     \
    para->index[5]                                                                      \
  );                                                                                    \
}                                                                                       \

#define KOKKOS_ATHREAD_REGISTER_REDUCE_1D(FUNC_NAME, ...)                               \
void FUNC_NAME##_1D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->reduce_result                                                                 \
  );                                                                                    \
}                                                                                       \
                                                                                       
#define KOKKOS_ATHREAD_REGISTER_REDUCE_2D(FUNC_NAME, ...)                               \
void FUNC_NAME##_2D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->reduce_result                                                                 \
  );                                                                                    \
}                                                                                       \
                                                                                       
#define KOKKOS_ATHREAD_REGISTER_REDUCE_3D(FUNC_NAME, ...)                               \
void FUNC_NAME##_3D(AthreadParamWrap* para) {                                           \
  const auto *f = reinterpret_cast<const __VA_ARGS__ *>(para->functor);                 \
  auto* F = const_cast<__VA_ARGS__ *>(f);                                               \
  F->operator() (                                                                       \
    para->index[0],                                                                     \
    para->index[1],                                                                     \
    para->index[2],                                                                     \
    para->reduce_result                                                                 \
  );                                                                                    \
}                                                                                       \

#endif // KOKKOS_ATHREAD_REGISTERMACROS_HPP_
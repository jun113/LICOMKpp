// __INCLUDE_FUNCTOR_HPP_START__
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_readyc.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_bclinc.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_icesnow.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_convadj.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_barotr.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_readyt.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_jra_daily.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_nextstep.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_tracer.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/kokkos-impl/kokkos_tracer.hpp"
#include "/home/export/base/shisuan/swiap/online/2024-SC/LICOMKpp-SC/src/util/pop_haloupdate.hpp"
// __INCLUDE_FUNCTOR_HPP_END__
#include "Kokkos_Athread_RegisterFunction.hpp"
#include "Kokkos_Athread_ParamWrap.h"

#include "simd.h"
#include "slave.h"

#include <string>

#include <time.h>

extern "C" void slave_cores_reduce(void* src_addr, void* dest_addr, 
		int units, int dtype, int optype, void *buf, int buf_item);


//ATHREAD DMA transfer definition 
__thread_local unsigned int D_COUNT  = 0;
__thread_local crts_rply_t  dma_rply = 0;
__thread_local crts_rply_t  l_rply   = 0;
__thread_local crts_rply_t  r_rply   = 0;

__thread_local_fix AthreadParamWrap kokkos_athread_local_param __attribute__ ((aligned(64)));

// Buffer for the operator of a reduction between slave cores
__thread_local double buf_reduce __attribute__ ((aligned(64)));

inline static int* arr_char_to_arr_int(const char *str, int &num_intv16) {
  // This function used SIMD for str cmp
  int key_len = 0;
  while (str[key_len] != '\0') {
    key_len += 1;
  }
  if (key_len > 32) {
    printf("len=%d, %s\n", key_len, str);
    printf("Current version the length of Functor's name must be less than 32.\n");
    exit(EXIT_FAILURE);
  }

  int *arr_int __attribute__ ((aligned(64))) = nullptr;

  if (key_len > 16) {
    arr_int = (int*)libc_aligned_malloc(sizeof(int) * 32);
  } else {
    arr_int = (int*)libc_aligned_malloc(sizeof(int) * 16);
  }

  for (int i = 0; i < key_len; ++i) {
    arr_int[i] = str[i];
  }

  if (key_len > 16) {
    for (int i = key_len; i < 32; ++i) {
      arr_int[i] = '0';
    }
    num_intv16 = 2;
  } else {
    for (int i = key_len; i < 16; ++i) {
      arr_int[i] = '0';
    }
    num_intv16 = 1;
  }
  return arr_int;
}

// TODO SIMD
inline static bool str_cmp(char* const str1, char* const str2) {
  int i1 = 0;
  int i2 = 0;
  while (str1[i1] != '\0' && str2[i2] != '\0') {
    if (str1[i1] != str2[i2]) {
      return true;
    }
    ++i1;
    ++i2;
  }
  return false;
}
extern int* g_athread_functor_key __attribute__ ((aligned(64)));

inline static bool str_cmp_intv16 (const intv16 &str1_simd, const intv16 &str2_simd, intv16 &result_simd) {

  // about 1 cycle
  result_simd = simd_vcmpeqw(str1_simd, str2_simd);

  // about 29 cycles
  if (simd_reduc_plusw(result_simd) == 16) {
    return true;
  } else {
    return false;
  }
}

extern AthradRegisterFunctionListNode* reg_func_list_node;
typedef void (*FunctionPointer)(AthreadParamWrap *);

static FunctionPointer lookup_fp () {
  using node = AthradRegisterFunctionListNode;

  intv16 *result_simd    = (intv16*)ldm_malloc(sizeof(intv16));
  intv16 *node_key_simd  = (intv16*)ldm_malloc(sizeof(intv16));
  intv16 *curr_key1_simd = (intv16*)ldm_malloc(sizeof(intv16));

  node* target_node = reg_func_list_node;

  if (kokkos_athread_local_param.num_intv16 == 1) {

    simd_load(curr_key1_simd[0], g_athread_functor_key);

    while (target_node != nullptr) {
      if (target_node->num_intv16 != kokkos_athread_local_param.num_intv16) {
        target_node = target_node->next;
      }
      simd_load(node_key_simd[0], target_node->key);
      if (str_cmp_intv16(node_key_simd[0], curr_key1_simd[0], result_simd[0])) {
        break;
      } else {
        target_node = target_node->next;
      }
    }
  } else if (kokkos_athread_local_param.num_intv16 == 2) {

    intv16 *curr_key2_simd = (intv16*)ldm_malloc(sizeof(intv16));

    simd_load(curr_key1_simd[0],   g_athread_functor_key);
    simd_load(curr_key2_simd[0], &(g_athread_functor_key[16]));

    while (target_node != nullptr) {

      if (target_node->num_intv16 != kokkos_athread_local_param.num_intv16) {
        target_node = target_node->next;
      }

      simd_load(node_key_simd[0], target_node->key);

      if (!str_cmp_intv16(node_key_simd[0], curr_key1_simd[0], result_simd[0])) {
        // compare first intv16
        target_node = target_node->next;
        continue;
      }

      simd_load(node_key_simd[0], &(target_node->key[16]));

      if (str_cmp_intv16(node_key_simd[0], curr_key2_simd[0], result_simd[0])) {
        // compare second intv16
        break;
      } else {
        target_node = target_node->next;
      }
    } // End while
    ldm_free(curr_key2_simd, sizeof(intv16));
  } else {
    printf ("error in %s, num_intv16 = %d\n", __PRETTY_FUNCTION__, kokkos_athread_local_param.num_intv16);
    exit(EXIT_FAILURE);
  }

  if (target_node == nullptr) {
    if (athread_tid == 0) {
      printf ("Error: \"");
      for (int i = 0; i < kokkos_athread_local_param.num_intv16; ++i) {
        for (int j = 0; j < 16; ++j) {
          printf("%c", g_athread_functor_key[i * 16 + j]);
        }
      }
      printf("\" was not found.\n");
    }
    exit(EXIT_FAILURE);
  }

  ldm_free(result_simd,    sizeof(intv16));
  ldm_free(node_key_simd,  sizeof(intv16));
  ldm_free(curr_key1_simd, sizeof(intv16));

  return target_node->fp;
}

extern "C" void register_kernel() {
  // __REGISTER_START__
  using node = AthradRegisterFunctionListNode;
  reg_func_list_node = (node*)libc_aligned_malloc(sizeof(node));
  node* curr_node __attribute__ ((aligned(64))) = reg_func_list_node;

  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc12", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc23", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc33", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc3_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc43", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc4_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc52", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc5_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc513", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc51_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc522", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc52_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc533", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc53_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc542", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc54_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc553", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc55_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc63", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc6_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc72", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc7_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc83", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc8_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc93", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc9_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc102", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc10_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc113", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc11_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvMomCen13", curr_node->num_intv16);
  curr_node->fp   = FuncAdvMomCen1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvMomFlu13", curr_node->num_intv16);
  curr_node->fp   = FuncAdvMomFlu1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvMomCen23", curr_node->num_intv16);
  curr_node->fp   = FuncAdvMomCen2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvMomFlu23", curr_node->num_intv16);
  curr_node->fp   = FuncAdvMomFlu2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc143", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc14_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc153", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc15_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc163", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc16_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc173", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc17_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc192", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc19_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc202", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc20_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc212", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc21_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc222", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc22_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc233", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc23_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc243", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc24_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyc253", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyc25_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc12", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc22", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc2_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc33", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc3_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc42", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc4_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc53", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc63", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc6_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc73", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc7_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc83", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc8_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc93", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc9_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc102", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc10_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc112", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc11_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc122", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc12_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc133", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc13_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc142", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc14_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc152", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc15_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc163", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc16_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc172", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc17_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc182", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc18_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBclinc193", curr_node->num_intv16);
  curr_node->fp   = FunctorBclinc19_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorIcesnow12", curr_node->num_intv16);
  curr_node->fp   = FunctorIcesnow1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorIcesnow23", curr_node->num_intv16);
  curr_node->fp   = FunctorIcesnow2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorConvadj12", curr_node->num_intv16);
  curr_node->fp   = FunctorConvadj1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorConvadj24", curr_node->num_intv16);
  curr_node->fp   = FunctorConvadj2_4D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorConvadj32", curr_node->num_intv16);
  curr_node->fp   = FunctorConvadj3_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr13", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr22", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr2_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr32", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr3_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr42", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr4_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr52", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr5_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr62", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr6_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr72", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr7_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr82", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr8_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr92", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr9_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr112", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr11_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr122", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr12_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr132", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr13_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr142", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr14_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr152", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr15_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr162", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr16_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr172", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr17_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr182", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr18_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr192", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr19_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr202", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr20_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorBarotr212", curr_node->num_intv16);
  curr_node->fp   = FunctorBarotr21_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt13", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt23", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt33", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt3_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt42", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt4_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt53", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt63", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt6_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt73", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt7_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt82", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt8_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt92", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt9_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt102", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt10_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt112", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt11_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt122", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt12_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorReadyt132", curr_node->num_intv16);
  curr_node->fp   = FunctorReadyt13_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorInterplationNearest2", curr_node->num_intv16);
  curr_node->fp   = FunctorInterplationNearest_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorNear12", curr_node->num_intv16);
  curr_node->fp   = FunctorNear1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorNear21", curr_node->num_intv16);
  curr_node->fp   = FunctorNear2_1D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorJRADaily12", curr_node->num_intv16);
  curr_node->fp   = FunctorJRADaily1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorJRADaily22", curr_node->num_intv16);
  curr_node->fp   = FunctorJRADaily2_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorJRADaily32", curr_node->num_intv16);
  curr_node->fp   = FunctorJRADaily3_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorNcarOceanFluxesJra2", curr_node->num_intv16);
  curr_node->fp   = FunctorNcarOceanFluxesJra_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorNextStep13", curr_node->num_intv16);
  curr_node->fp   = FunctorNextStep1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorNextStep22", curr_node->num_intv16);
  curr_node->fp   = FunctorNextStep2_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer12", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer23", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer32", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer3_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer43", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer4_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer53", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer62", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer6_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer73", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer7_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer83", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer8_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraCenTsp13", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraCenTsp1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp23", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp33", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp3_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp42", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp4_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp53", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp62", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp6_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer153", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer15_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer163", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer16_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer173", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer17_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer183", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer18_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer192", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer19_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer203", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer20_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer212", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer21_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer222", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer22_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer233", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer23_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer242", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer24_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer252", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer25_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer262", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer26_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer272", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer27_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer282", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer28_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer292", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer29_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer303", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer30_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer313", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer31_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer322", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer32_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer332", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer33_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer343", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer34_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer353", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer35_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer12", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer1_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer23", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer32", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer3_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer43", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer4_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer53", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer62", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer6_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer73", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer7_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer83", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer8_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraCenTsp13", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraCenTsp1_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp23", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp2_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp33", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp3_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp42", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp4_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp53", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp5_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FuncAdvTraTsp62", curr_node->num_intv16);
  curr_node->fp   = FuncAdvTraTsp6_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer153", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer15_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer163", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer16_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer173", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer17_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer183", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer18_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer192", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer19_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer203", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer20_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer212", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer21_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer222", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer22_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer233", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer23_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer242", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer24_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer252", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer25_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer262", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer26_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer272", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer27_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer282", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer28_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer292", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer29_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer303", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer30_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer313", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer31_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer322", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer32_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer332", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer33_2D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer343", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer34_3D;
  curr_node->next = (node*)libc_aligned_malloc(sizeof(node));

  curr_node       = curr_node->next;
  curr_node->key  = arr_char_to_arr_int (
      "FunctorTracer353", curr_node->num_intv16);
  curr_node->fp   = FunctorTracer35_3D;
  curr_node->next = nullptr;
  // __REGISTER_END__
  return ;
}

extern "C" void parallel_for_1D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start = kokkos_athread_local_param.range[0][0];
  const int end   = kokkos_athread_local_param.range[0][1];
  const int len   = end - start;
  const int times = (len + ATHREAD_SLAVE_CORES - 1) / ATHREAD_SLAVE_CORES;

  for (int i = 0; i < times; ++i) {
    int index = start + athread_tid * times + i;
    if (index < end) {
      kokkos_athread_local_param.index[0] = index;
      fp(&kokkos_athread_local_param);
    }
  }
  return ;
}

extern "C" void parallel_for_2D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start0 = kokkos_athread_local_param.range[0][0];
  const int end0   = kokkos_athread_local_param.range[0][1];
  const int start1 = kokkos_athread_local_param.range[1][0];
  const int end1   = kokkos_athread_local_param.range[1][1];
  const int tile0  = kokkos_athread_local_param.tile[0];
  const int tile1  = kokkos_athread_local_param.tile[1];

  const int num_tiles0 = (end0 - start0 + tile0 - 1) / tile0;
  const int num_tiles1 = (end1 - start1 + tile1 - 1) / tile1;

  const int num_tiles = num_tiles0 * num_tiles1;
  
  for (int index_tile = athread_tid; index_tile < num_tiles; 
      index_tile += ATHREAD_SLAVE_CORES) {
    const int index_tile0 = index_tile / num_tiles1;
    const int index_tile1 = index_tile % num_tiles1;

    for (int index00 = 0; index00 < tile0; ++ index00) {
      const int index0 = start0 + index_tile0 * tile0 + index00;
      if (index0 >= end0) { break; }
      for (int index_11 = 0; index_11 < tile1; ++ index_11) {
        const int index1 = start1 + index_tile1 * tile1 + index_11;
        if (index1 >= end1) { break; }
        kokkos_athread_local_param.index[0] = index0;
        kokkos_athread_local_param.index[1] = index1;
        fp(&kokkos_athread_local_param);
      }
    }
  }
  return ;
}

extern "C" void parallel_for_3D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start[3] = {kokkos_athread_local_param.range[0][0], kokkos_athread_local_param.range[1][0], kokkos_athread_local_param.range[2][0]};
  const int end[3]   = {kokkos_athread_local_param.range[0][1], kokkos_athread_local_param.range[1][1], kokkos_athread_local_param.range[2][1]};
  const int tile[3]  = {kokkos_athread_local_param.tile[0],     kokkos_athread_local_param.tile[1],     kokkos_athread_local_param.tile[2]};

  const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];
  const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];
  const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];

  const int tmp0 = num_tiles1 * num_tiles2;
  const int num_tiles = num_tiles0 * tmp0;
  
  for (int index_tile = athread_tid; index_tile < num_tiles; index_tile += ATHREAD_SLAVE_CORES) {

    const int index_tile0 = index_tile / tmp0;
    const int tmp1        = index_tile % tmp0;
    const int index_tile1 = tmp1      / num_tiles2;
    const int index_tile2 = tmp1      % num_tiles2;

    for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {
      const int index0 = start[0] + index_tile0 * tile[0] + index_00;
      if (index0 >= end[0]) { break; }
      for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {
        const int index1 = start[1] + index_tile1 * tile[1] + index_11;
        if (index1 >= end[1]) { break; }
        for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {
          const int index2 = start[2] + index_tile2 * tile[2] + index_22;
          if (index2 >= end[2]) { break; }

          kokkos_athread_local_param.index[0] = index0;
          kokkos_athread_local_param.index[1] = index1;
          kokkos_athread_local_param.index[2] = index2;
          fp(&kokkos_athread_local_param);
        }
      }
    }
  }
  return ;
}

extern "C" void parallel_for_4D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start[4] = {kokkos_athread_local_param.range[0][0], kokkos_athread_local_param.range[1][0], 
      kokkos_athread_local_param.range[2][0], kokkos_athread_local_param.range[3][0]};
  const int end[4]   = {kokkos_athread_local_param.range[0][1], kokkos_athread_local_param.range[1][1], 
      kokkos_athread_local_param.range[2][1], kokkos_athread_local_param.range[3][1]};
  const int tile[4]  = {kokkos_athread_local_param.tile[0],     kokkos_athread_local_param.tile[1],     
      kokkos_athread_local_param.tile[2],     kokkos_athread_local_param.tile[3]};

  const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];
  const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];
  const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];
  const int num_tiles3 = (end[3] - start[3] + tile[3] - 1) / tile[3];

  const int tmp0      = num_tiles2 * num_tiles3;
  const int tmp1      = num_tiles1 * tmp0;
  const int num_tiles = num_tiles0 * tmp1;
  
  for (int index_tile = athread_tid; index_tile < num_tiles; index_tile += ATHREAD_SLAVE_CORES) {

    const int index_tile0 = index_tile / tmp1;
    const int tmp2        = index_tile % tmp1;
    const int index_tile1 = tmp2       / tmp0;
    const int tmp3        = tmp2       % tmp0;
    const int index_tile2 = tmp3       / num_tiles3;
    const int index_tile3 = tmp3       % num_tiles3;

    for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {
      const int index0 = start[0] + index_tile0 * tile[0] + index_00;
      if (index0 >= end[0]) { break; }
      for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {
        const int index1 = start[1] + index_tile1 * tile[1] + index_11;
        if (index1 >= end[1]) { break; }
        for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {
          const int index2 = start[2] + index_tile2 * tile[2] + index_22;
          if (index2 >= end[2]) { break; }
          for (int index_33 = 0; index_33 < tile[3]; ++ index_33) {
            const int index3 = start[3] + index_tile3 * tile[3] + index_33;
            if (index3 >= end[3]) { break; }

            kokkos_athread_local_param.index[0] = index0;
            kokkos_athread_local_param.index[1] = index1;
            kokkos_athread_local_param.index[2] = index2;
            kokkos_athread_local_param.index[3] = index3;
            fp(&kokkos_athread_local_param);
          }
        }
      }
    }
  }
  return ;
}

extern "C" void parallel_reduce_1D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start = kokkos_athread_local_param.range[0][0];
  const int end   = kokkos_athread_local_param.range[0][1];
  const int len   = end - start;
  const int times = (len + ATHREAD_SLAVE_CORES - 1) / ATHREAD_SLAVE_CORES;

  for (int i = 0; i < times; ++i) {
    int index = start + athread_tid * times + i;
    if (index < end) {
      kokkos_athread_local_param.index[0] = index;
      fp(&kokkos_athread_local_param);
    }
  }
  // reduce from all slave cores
  slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
      &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
  if (athread_tid == 0) {
    host_param->reduce_result = kokkos_athread_local_param.reduce_result;
  }
  return ;
}

extern "C" void parallel_reduce_2D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start0 = kokkos_athread_local_param.range[0][0];
  const int end0   = kokkos_athread_local_param.range[0][1];
  const int start1 = kokkos_athread_local_param.range[1][0];
  const int end1   = kokkos_athread_local_param.range[1][1];
  const int tile0  = kokkos_athread_local_param.tile[0];
  const int tile1  = kokkos_athread_local_param.tile[1];

  const int num_tiles0 = (end0 - start0 + tile0 - 1) / tile0;
  const int num_tiles1 = (end1 - start1 + tile1 - 1) / tile1;

  const int num_tiles = num_tiles0 * num_tiles1;
  
  for (int index_tile = athread_tid; index_tile < num_tiles; 
      index_tile += ATHREAD_SLAVE_CORES) {
    const int index_tile0 = index_tile / num_tiles1;
    const int index_tile1 = index_tile % num_tiles1;

    for (int index00 = 0; index00 < tile0; ++ index00) {
      const int index0 = start0 + index_tile0 * tile0 + index00;
      if (index0 >= end0) { break; }
      for (int index_11 = 0; index_11 < tile1; ++ index_11) {
        const int index1 = start1 + index_tile1 * tile1 + index_11;
        if (index1 >= end1) { break; }
        kokkos_athread_local_param.index[0] = index0;
        kokkos_athread_local_param.index[1] = index1;
        fp(&kokkos_athread_local_param);
      }
    }
  }
  // reduce from all slave cores
  slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
      &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
  if (athread_tid == 0) {
    host_param->reduce_result = kokkos_athread_local_param.reduce_result;
  }
  return ;
}

extern "C" void parallel_reduce_3D(AthreadParamWrap* host_param) {

  athread_dma_iget(&kokkos_athread_local_param, host_param, sizeof(AthreadParamWrap), &dma_rply);
  D_COUNT++;
  athread_dma_wait_value(&dma_rply, D_COUNT);

  const FunctionPointer fp __attribute__ ((aligned(64))) = lookup_fp();

  const int ATHREAD_SLAVE_CORES = 64;

  const int start[3] = {kokkos_athread_local_param.range[0][0], kokkos_athread_local_param.range[1][0], kokkos_athread_local_param.range[2][0]};
  const int end[3]   = {kokkos_athread_local_param.range[0][1], kokkos_athread_local_param.range[1][1], kokkos_athread_local_param.range[2][1]};
  const int tile[3]  = {kokkos_athread_local_param.tile[0],     kokkos_athread_local_param.tile[1],     kokkos_athread_local_param.tile[2]};

  const int num_tiles0 = (end[0] - start[0] + tile[0] - 1) / tile[0];
  const int num_tiles1 = (end[1] - start[1] + tile[1] - 1) / tile[1];
  const int num_tiles2 = (end[2] - start[2] + tile[2] - 1) / tile[2];

  const int tmp0 = num_tiles1 * num_tiles2;
  const int num_tiles = num_tiles0 * tmp0;
  
  for (int index_tile = athread_tid; index_tile < num_tiles; index_tile += ATHREAD_SLAVE_CORES) {

    const int index_tile0 = index_tile / tmp0;
    const int tmp1        = index_tile % tmp0;
    const int index_tile1 = tmp1      / num_tiles2;
    const int index_tile2 = tmp1      % num_tiles2;

    for (int index_00 = 0; index_00 < tile[0]; ++ index_00) {
      const int index0 = start[0] + index_tile0 * tile[0] + index_00;
      if (index0 >= end[0]) { break; }
      for (int index_11 = 0; index_11 < tile[1]; ++ index_11) {
        const int index1 = start[1] + index_tile1 * tile[1] + index_11;
        if (index1 >= end[1]) { break; }
        for (int index_22 = 0; index_22 < tile[2]; ++ index_22) {
          const int index2 = start[2] + index_tile2 * tile[2] + index_22;
          if (index2 >= end[2]) { break; }

          kokkos_athread_local_param.index[0] = index0;
          kokkos_athread_local_param.index[1] = index1;
          kokkos_athread_local_param.index[2] = index2;
          fp(&kokkos_athread_local_param);
        }
      }
    }
  }
  // reduce from all slave cores
  slave_cores_reduce(&(kokkos_athread_local_param.reduce_result),
      &(kokkos_athread_local_param.reduce_result), 1, athread_double, OP_add, &buf_reduce, 1);
  if (athread_tid == 0) {
    host_param->reduce_result = kokkos_athread_local_param.reduce_result;
  }
  return ;
}
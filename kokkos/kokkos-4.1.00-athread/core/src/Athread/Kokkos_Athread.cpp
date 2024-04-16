//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#if defined(__sw_host__)
#include <Kokkos_Core.hpp>
#endif // defined(__sw_host__)

#include <Athread/Kokkos_Athread.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

#include "athread.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

#include <string>
/*--------------------------------------------------------------------------*/
extern "C" void launch_register_kernel();

// TODO memory free
AthradRegisterFunctionListNode* reg_func_list_node __attribute__ ((aligned(64))) = nullptr;
// DMA: 128B, sizeof(uint) = 4B
// 128 / 4 = 32
int* g_athread_functor_key __attribute__ ((aligned(64))) = 
    (int *)libc_aligned_malloc(32 * sizeof(int));
// uintv16: 16 * uint
// uintv16* g_athread_functor_key_simd __attribute__ ((aligned(64))) = 
//     (uintv16 *)libc_aligned_malloc(2 * sizeof(uintv16));

namespace Kokkos {
namespace Impl {

bool AthreadInternal::is_initialized() { return m_is_initialized; }

void AthreadInternal::initialize() {
  if (is_initialized()) return;

  athread_init();

  // TODO
  // get_current_max_cores();

  launch_register_kernel();

  Impl::SharedAllocationRecord<void, void>::tracking_enable();

  m_is_initialized = true;
}

void AthreadInternal::finalize() {
  if (m_thread_team_data.scratch_buffer()) {
    m_thread_team_data.disband_team();
    m_thread_team_data.disband_pool();

    Kokkos::HostSpace space;

    space.deallocate(m_thread_team_data.scratch_buffer(),
                     m_thread_team_data.scratch_bytes());

    m_thread_team_data.scratch_assign(nullptr, 0, 0, 0, 0, 0);
  }

  Kokkos::Profiling::finalize();

  athread_halt();
  m_is_initialized = false;
}

AthreadInternal& AthreadInternal::singleton() {
  static AthreadInternal* self = nullptr;
  if (!self) {
    self = new AthreadInternal();
  }
  return *self;
}

// Resize thread team data scratch memory
void AthreadInternal::resize_thread_team_data(size_t pool_reduce_bytes,
                                             size_t team_reduce_bytes,
                                             size_t team_shared_bytes,
                                             size_t thread_local_bytes) {
  if (pool_reduce_bytes < 512) pool_reduce_bytes = 512;
  if (team_reduce_bytes < 512) team_reduce_bytes = 512;

  const size_t old_pool_reduce  = m_thread_team_data.pool_reduce_bytes();
  const size_t old_team_reduce  = m_thread_team_data.team_reduce_bytes();
  const size_t old_team_shared  = m_thread_team_data.team_shared_bytes();
  const size_t old_thread_local = m_thread_team_data.thread_local_bytes();
  const size_t old_alloc_bytes  = m_thread_team_data.scratch_bytes();

  // Allocate if any of the old allocation is tool small:

  const bool allocate = (old_pool_reduce < pool_reduce_bytes) ||
                        (old_team_reduce < team_reduce_bytes) ||
                        (old_team_shared < team_shared_bytes) ||
                        (old_thread_local < thread_local_bytes);

  if (allocate) {
    Kokkos::HostSpace space;

    if (old_alloc_bytes) {
      m_thread_team_data.disband_team();
      m_thread_team_data.disband_pool();

      space.deallocate("Kokkos::Athread::scratch_mem",
                       m_thread_team_data.scratch_buffer(),
                       m_thread_team_data.scratch_bytes());
    }

    if (pool_reduce_bytes < old_pool_reduce) {
      pool_reduce_bytes = old_pool_reduce;
    }
    if (team_reduce_bytes < old_team_reduce) {
      team_reduce_bytes = old_team_reduce;
    }
    if (team_shared_bytes < old_team_shared) {
      team_shared_bytes = old_team_shared;
    }
    if (thread_local_bytes < old_thread_local) {
      thread_local_bytes = old_thread_local;
    }

    const size_t alloc_bytes =
        HostThreadTeamData::scratch_size(pool_reduce_bytes, team_reduce_bytes,
                                         team_shared_bytes, thread_local_bytes);

    void* ptr = nullptr;
    try {
      ptr = space.allocate("Kokkos::Athread::scratch_mem", alloc_bytes);
    } catch (Kokkos::Experimental::RawMemoryAllocationFailure const& failure) {
      // For now, just rethrow the error message the existing way
      Kokkos::Impl::throw_runtime_exception(failure.get_error_message());
    }

    m_thread_team_data.scratch_assign(static_cast<char*>(ptr), alloc_bytes,
                                      pool_reduce_bytes, team_reduce_bytes,
                                      team_shared_bytes, thread_local_bytes);

    HostThreadTeamData* pool[1] = {&m_thread_team_data};

    m_thread_team_data.organize_pool(pool, 1);
    m_thread_team_data.organize_team(1);
  }
}

int AthreadInternal::get_current_max_cores() {
  unsigned int max_cores = 0;
  athread_res_t* res_curr = new athread_res_t;

  if (res_curr) {
    athread_res_getcurrent(res_curr);
    int err = 0;
    err = athread_res_getpenum(res_curr, &max_cores);
    delete res_curr;
    if (err != 0) {
      // TODO Athread error
      return 0;
    }
  }
  res_curr = nullptr;
  g_athread_hardware_max_cores = max_cores;
  return max_cores;
}
}  // namespace Impl

Athread::Athread()
    : m_space_instance(&Impl::AthreadInternal::singleton(),
                       [](Impl::AthreadInternal*) {}) {}

void Athread::print_configuration(std::ostream& os, bool /*verbose*/) const {
  os << "Host Athread Execution Space:\n";
  os << "  KOKKOS_ENABLE_ATHREAD: yes\n";

#ifdef KOKKOS_INTERNAL_NOT_PARALLEL
  os << "Kokkos atomics disabled\n";
#endif

  os << "\nAthread Runtime Configuration:\n";

  printf("Host Athread Execution Space:\n");
  printf("KOKKOS_ENABLE_ATHREAD: yes\n");
}

bool Athread::impl_is_initialized() {
  return Impl::AthreadInternal::singleton().is_initialized();
}

void Athread::impl_initialize(InitializationSettings const&) {
  Impl::AthreadInternal::singleton().initialize();
}

void Athread::impl_finalize() { Impl::AthreadInternal::singleton().finalize(); }

const char* Athread::name() { return "Athread"; }

namespace Impl {

#if defined(__sw_host__)
int g_serial_space_factory_initialized =
    initialize_space_factory<Athread>("100_Athread");
#endif // defined(__sw_host__)


}  // namespace Impl

}  // namespace Kokkos

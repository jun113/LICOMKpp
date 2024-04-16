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

#ifndef KOKKO_ATHREAD_PARALLEL_MDRANGE_HPP
#define KOKKO_ATHREAD_PARALLEL_MDRANGE_HPP

#include <athread.h>
#include <Kokkos_Parallel.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

#include <iostream>

namespace Kokkos {
namespace Impl {

// template <class FunctorType, class... Traits>
// class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
//                   Kokkos::Athread> {
//  private:
//   using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
//   using Policy        = typename MDRangePolicy::impl_range_policy;

//   using iterate_type = typename Kokkos::Impl::HostIterateTile<
//       MDRangePolicy, FunctorType, typename MDRangePolicy::work_tag, void>;

//   const iterate_type m_iter;

//   void exec() const {
//     const typename Policy::member_type e = m_iter.m_rp.m_num_tiles;
//     for (typename Policy::member_type i = 0; i < e; ++i) {
//       m_iter(i);
//     }
//   }

//  public:
//   inline void execute() const { this->exec(); }
//   template <typename Policy, typename Functor>
//   static int max_tile_size_product(const Policy&, const Functor&) {
//     /**
//      * 1024 here is just our guess for a reasonable max tile size,
//      * it isn't a hardware constraint. If people see a use for larger
//      * tile size products, we're happy to change this.
//      */
//     return 1024;
//   }
//   inline ParallelFor(const FunctorType& arg_functor,
//                      const MDRangePolicy& arg_policy)
//       : m_iter(arg_policy, arg_functor) {}
// };

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Athread> {
 private:
  using MDPolicy = Kokkos::MDRangePolicy<Traits...>;

  using array_index_type = typename MDPolicy::array_index_type;
  using index_type       = typename MDPolicy::index_type;
  using LaunchBounds     = typename MDPolicy::launch_bounds;

  const FunctorType m_functor;
  const MDPolicy    m_policy;

  void exec() const {
    if (m_policy.m_num_tiles == 0) return;

    static AthreadParamWrap param;

    if (MDPolicy::rank == 2 ) {

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];

      athread_get_key<FunctorType>(2, param.num_intv16);
      param.functor = reinterpret_cast<const void *>(&m_functor);
      athread_parallel_for_launch_2D(&param);
    } else if (MDPolicy::rank == 3) {

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];
      // printf ("0: %d %d, 1: %d %d, 2: %d %d\n", 
      //     m_policy.m_lower[0], m_policy.m_upper[0],
      //     m_policy.m_lower[1], m_policy.m_upper[1],
      //     m_policy.m_lower[2], m_policy.m_upper[2]);

      athread_get_key<FunctorType>(3, param.num_intv16);
      param.functor = reinterpret_cast<const void *>(&m_functor);
      athread_parallel_for_launch_3D(&param);
    } else if (MDPolicy::rank == 4) {

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.range[3][0] = m_policy.m_lower[3];
      param.range[3][1] = m_policy.m_upper[3];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];
      param.tile[3]     = m_policy.m_tile[3];

      athread_get_key<FunctorType>(4, param.num_intv16);
      param.functor = reinterpret_cast<const void *>(&m_functor);
      athread_parallel_for_launch_4D(&param);
    } else if (MDPolicy::rank == 5) {
      std::cout<<"TODO MDRangePolicy::rank == 5"<<std::endl;
      return ;
    } else if (MDPolicy::rank == 6) {
      std::cout<<"TODO MDRangePolicy::rank == 6"<<std::endl;
      return ;
    }
  }

 public:
  inline void execute() const { this->exec(); }
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const MDPolicy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline ParallelFor(const FunctorType& arg_functor,
                     const MDPolicy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::MDRangePolicy<Traits...>, Kokkos::Athread> {
 private:
  using MDRangePolicy = Kokkos::MDRangePolicy<Traits...>;
  using Policy        = typename MDRangePolicy::impl_range_policy;
  using FunctorType   = typename CombinedFunctorReducerType::functor_type;
  using ReducerType   = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename MDRangePolicy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using value_type     = typename ReducerType::value_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const MDRangePolicy m_policy;
  const pointer_type m_result_ptr;

  inline void exec(reference_type update) const {

    if (m_policy.m_num_tiles == 0) return;

    static AthreadParamWrap param;

    if (MDRangePolicy::rank == 2 ) {

      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];

      athread_get_key<FunctorType>(2, param.num_intv16);

      param.reduce_result = static_cast<double>(update);
      param.functor = reinterpret_cast<const void *>(
          &(m_functor_reducer.get_functor()));

      athread_parallel_reduce_launch_2D(&param);

    } else if (MDRangePolicy::rank == 3) {
      param.range[0][0] = m_policy.m_lower[0];
      param.range[0][1] = m_policy.m_upper[0];
      param.range[1][0] = m_policy.m_lower[1];
      param.range[1][1] = m_policy.m_upper[1];
      param.range[2][0] = m_policy.m_lower[2];
      param.range[2][1] = m_policy.m_upper[2];
      param.tile[0]     = m_policy.m_tile[0];
      param.tile[1]     = m_policy.m_tile[1];
      param.tile[2]     = m_policy.m_tile[2];

      athread_get_key<FunctorType>(3, param.num_intv16);
      param.reduce_result = static_cast<double>(update);

      param.functor = reinterpret_cast<const void *>(&(m_functor_reducer.get_functor()));
      athread_parallel_reduce_launch_3D(&param);
    } else if (MDRangePolicy::rank == 4) {
      std::cout<<"TODO MDRangePolicy::rank == 4"<<std::endl;
    } else if (MDRangePolicy::rank == 5) {
      std::cout<<"TODO MDRangePolicy::rank == 5"<<std::endl;
    } else if (MDRangePolicy::rank == 6) {
      std::cout<<"TODO MDRangePolicy::rank == 6"<<std::endl;
    }

    update = static_cast<reference_type>(param.reduce_result);
  }

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    /**
     * 1024 here is just our guess for a reasonable max tile size,
     * it isn't a hardware constraint. If people see a use for larger
     * tile size products, we're happy to change this.
     */
    return 1024;
  }
  inline void execute() const {
    const ReducerType& reducer = m_functor_reducer.get_reducer();
    if (m_result_ptr == nullptr) {
      printf("error in ParallelReduce, execute(), m_result_ptr == nullptr\n");
      exit(0);
    }
    pointer_type ptr = m_result_ptr;
    reference_type update = reducer.init(ptr);
    this->exec(update);
    reducer.final(ptr);
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const MDRangePolicy& arg_policy,
                 const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<ViewType>::value,
                  "Kokkos::Athread reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Athread reduce result must be a View accessible from "
        "HostSpace");
  }
};


}  // namespace Impl
}  // namespace Kokkos

#endif

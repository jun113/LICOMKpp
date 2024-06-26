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

#ifndef KOKKOS_ATHREAD_ZEROMEMSET_HPP
#define KOKKOS_ATHREAD_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>
#include <Athread/Kokkos_Athread.hpp>

#include <type_traits>

namespace Kokkos {
namespace Impl {

// We only need to provide a specialization for Athread if there is a host
// parallel execution space since the specialization for
// DefaultHostExecutionSpace is defined elsewhere.
struct DummyExecutionSpace;
template <class T, class... P>
struct ZeroMemset<
    std::conditional_t<!std::is_same<Athread, DefaultHostExecutionSpace>::value,
                       Athread, DummyExecutionSpace>,
    View<T, P...>>
    : public ZeroMemset<DefaultHostExecutionSpace, View<T, P...>> {
  using Base = ZeroMemset<DefaultHostExecutionSpace, View<T, P...>>;
  using Base::Base;

  ZeroMemset(const Athread&, const View<T, P...>& dst,
             typename View<T, P...>::const_value_type& value)
      : Base(dst, value) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_ATHREAD_ZEROMEMSET_HPP)

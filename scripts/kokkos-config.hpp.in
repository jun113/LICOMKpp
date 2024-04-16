#ifndef LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CONFIG_
#define LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CONFIG_

#include "def-undef.h"
#ifdef LICOM_ENABLE_KOKKOS

#include <Kokkos_Core.hpp>

//using Layout = Kokkos::LayoutLeft;
using Layout = Kokkos::LayoutRight;

using koArr2D = Kokkos::Array<int64_t, 2>;
using koArr3D = Kokkos::Array<int64_t, 3>;
using koArr4D = Kokkos::Array<int64_t, 4>;

#if (defined KOKKOS_ENABLE_CUDA) || (defined KOKKOS_ENABLE_HIP)
#define KOKKOS_ENABLE_DEVICE_MEM_SPACE
#else
#define KOKKOS_ENABLE_HOST_MEM_SPACE
#endif

#ifdef KOKKOS_ENABLE_DEVICE_MEM_SPACE
using ViewInt1D = Kokkos::View<int*,     Layout>;
using ViewInt2D = Kokkos::View<int**,    Layout>;
using ViewInt3D = Kokkos::View<int***,   Layout>;
using ViewInt4D = Kokkos::View<int****,  Layout>;
using ViewInt5D = Kokkos::View<int*****, Layout>;

using ViewFloat1D = Kokkos::View<float*,     Layout>;
using ViewFloat2D = Kokkos::View<float**,    Layout>;
using ViewFloat3D = Kokkos::View<float***,   Layout>;
using ViewFloat4D = Kokkos::View<float****,  Layout>;
using ViewFloat5D = Kokkos::View<float*****, Layout>;

using ViewDouble1D = Kokkos::View<double*,     Layout>;
using ViewDouble2D = Kokkos::View<double**,    Layout>;
using ViewDouble3D = Kokkos::View<double***,   Layout>;
using ViewDouble4D = Kokkos::View<double****,  Layout>;
using ViewDouble5D = Kokkos::View<double*****, Layout>;

using UnManagedViewInt1D = Kokkos::View<int*,     Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewInt2D = Kokkos::View<int**,    Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewInt3D = Kokkos::View<int***,   Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewInt4D = Kokkos::View<int****,  Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewInt5D = Kokkos::View<int*****, Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

using UnManagedViewFloat1D = Kokkos::View<float*,     Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewFloat2D = Kokkos::View<float**,    Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewFloat3D = Kokkos::View<float***,   Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewFloat4D = Kokkos::View<float****,  Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewFloat5D = Kokkos::View<float*****, Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

using UnManagedViewDouble1D = Kokkos::View<double*,     Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewDouble2D = Kokkos::View<double**,    Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewDouble3D = Kokkos::View<double***,   Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewDouble4D = Kokkos::View<double****,  Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
using UnManagedViewDouble5D = Kokkos::View<double*****, Layout,
		Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

#else

using ViewInt1D = Kokkos::View<int*,     Layout, Kokkos::MemoryUnmanaged>;
using ViewInt2D = Kokkos::View<int**,    Layout, Kokkos::MemoryUnmanaged>;
using ViewInt3D = Kokkos::View<int***,   Layout, Kokkos::MemoryUnmanaged>;
using ViewInt4D = Kokkos::View<int****,  Layout, Kokkos::MemoryUnmanaged>;
using ViewInt5D = Kokkos::View<int*****, Layout, Kokkos::MemoryUnmanaged>;

using ViewFloat1D = Kokkos::View<float*,     Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat2D = Kokkos::View<float**,    Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat3D = Kokkos::View<float***,   Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat4D = Kokkos::View<float****,  Layout, Kokkos::MemoryUnmanaged>;
using ViewFloat5D = Kokkos::View<float*****, Layout, Kokkos::MemoryUnmanaged>;

using ViewDouble1D = Kokkos::View<double*,     Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble2D = Kokkos::View<double**,    Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble3D = Kokkos::View<double***,   Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble4D = Kokkos::View<double****,  Layout, Kokkos::MemoryUnmanaged>;
using ViewDouble5D = Kokkos::View<double*****, Layout, Kokkos::MemoryUnmanaged>;
#endif

// tile

// static auto tile2D = koArr2D{4, 16};
// static auto tile3D = koArr3D{1, 4, 16};
// CUDA tile
// static auto tile2D = koArr2D{2, 16};
// static auto tile3D = koArr3D{2, 4, 32};

// static auto tile2D = koArr2D{2, 16};
// static auto tile3D = koArr3D{2, 4, 16};

// static auto tile2D = koArr2D{4, 16};
// static auto tile3D = koArr3D{2, 4, 32};
//
// HIP tile
// static auto tile2D = koArr2D{4, 16};
// static auto tile3D = koArr3D{2, 4, 32};
// 1.662 10.909
//static auto tile2D = koArr2D{2, 16};
//static auto tile3D = koArr3D{2, 4, 32};
// 1.721 9.751
//static auto tile2D = koArr2D{4, 16};
//static auto tile3D = koArr3D{2, 4, 16};
// 1.708 10.290
//static auto tile2D = koArr2D{4, 32};
//static auto tile3D = koArr3D{2, 4, 32};
//
// Sunway tile
//static auto tile2D = koArr2D{6, 16};
//static auto tile3D = koArr3D{5, 16, 32};
// SYPD 0.56
// static auto tile2D = koArr2D{1, 94};
// static auto tile3D = koArr3D{1, 1, 94};
// SYPD 0.68
// static auto tile2D = koArr2D{1, 25};
// static auto tile3D = koArr3D{1, 1, 25};

// static auto tile2D = koArr2D{1, 124};
// static auto tile3D = koArr3D{1, 1, 124};

// static auto tile2D = koArr2D{2, 62};
// static auto tile3D = koArr3D{1, 2, 62};

static auto tile2D = koArr2D{4, 16};
static auto tile3D = koArr3D{4, 4, 64};


#endif // LICOM_ENABLE_KOKKOS
#endif //LICOM3_KOKKOS_SRC_HEAD_KOKKOS_CONFIG_

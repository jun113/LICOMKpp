#ifndef KOKKOS_ATHREAD_PARAMWRAP_H_
#define KOKKOS_ATHREAD_PARAMWRAP_H_

struct AthreadParamWrap {
	int         range[6][2];
	int         tile[6];
	int         index[6];
	const void* functor;
	double      reduce_result;
	int         num_intv16;
};

// struct AthreadForParamWrap1D {
// 	int   range[1][2];
// 	int   index[1];
// 	char* key;
// 	const void* functor;
// };

// struct AthreadForParamWrap2D {
// 	int   range[2][2];
// 	int   tile[2];
// 	int   index[2];
// 	char* key;
// 	const void* functor;
// };

// struct AthreadForParamWrap3D {
// 	int   range[3][2];
// 	int   tile[3];
// 	int   index[3];
// 	char* key;
// 	const void* functor;
// };

// struct AthreadForParamWrap4D {
// 	int   range[4][2];
// 	int   tile[4];
// 	int   index[4];
// 	char* key;
// 	const void* functor;
// };

// struct AthreadForParamWrap5D {
// 	int   range[5][2];
// 	int   tile[5];
// 	int   index[5];
// 	char* key;
// 	const void* functor;
// };

// struct AthreadForParamWrap6D {
// 	int   range[6][2];
// 	int   tile[6];
// 	int   index[6];
// 	char* key;
// 	const void* functor;
// };

struct AthreadForParamWrap1D {
	int   range[1][2];
	char* key;
	const void* functor;
};

struct AthreadForParamWrap2D {
	int   range[2][2];
	int   tile[2];
	char* key;
	const void* functor;
};

struct AthreadForParamWrap3D {
	int   range[3][2];
	int   tile[3];
	char* key;
	const void* functor;
};

struct AthreadForParamWrap4D {
	int   range[4][2];
	int   tile[4];
	char* key;
	const void* functor;
};

struct AthreadForParamWrap5D {
	int   range[5][2];
	int   tile[5];
	char* key;
	const void* functor;
};

struct AthreadForParamWrap6D {
	int   range[6][2];
	int   tile[6];
	char* key;
	const void* functor;
};

struct AthreadReduceParamWrap1D {
	int   range[1][2];
	double reduce_result;
	char* key;
	const void* functor;
};

struct AthreadReduceParamWrap2D {
	int   range[2][2];
	int   tile[2];
	double* reduce_result;
	char*   key;
	const void* functor;
};

struct AthreadReduceParamWrap3D {
	int   range[3][2];
	int   tile[3];
	double reduce_result;
	char* key;
	const void* functor;
};

struct AthreadReduceParamWrap4D {
	int   range[4][2];
	int   tile[4];
	double reduce_result;
	char* key;
	const void* functor;
};

struct AthreadReduceParamWrap5D {
	int   range[5][2];
	int   tile[5];
	double reduce_result;
	char* key;
	const void* functor;
};

struct AthreadReduceParamWrap6D {
	int   range[6][2];
	int   tile[6];
	double reduce_result;
	char* key;
	const void* functor;
};

struct AthreadLaunchParamWrap {
	int    index[6];
	double reduce_result;
	const void* functor;
};



#endif // KOKKOS_ATHREAD_PARAMWRAP_H_
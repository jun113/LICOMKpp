#ifdef __sw_slave__

#include "athread_param.h"

#include "slave.h"

//ATHREAD DMA transfer definition 
// __thread_local unsigned int D_COUNT  = 0;
// __thread_local crts_rply_t  dma_rply = 0;
// __thread_local crts_rply_t  l_rply   = 0;
// __thread_local crts_rply_t  r_rply   = 0;

__thread_local struct HaloTransposeDouble3D local_param __attribute__ ((aligned(64)));

void athread_get_halo_transpose_double_slave (struct HaloTransposeDouble3D* arg_param) {

	struct HaloTransposeDouble3D local_param;
  athread_dma_get(&local_param, arg_param, 
      sizeof(struct HaloTransposeDouble3D));
  // athread_dma_iget(&local_param, arg_param, sizeof(struct HaloTransposeDouble3D), &dma_rply);
  // D_COUNT++;
  // athread_dma_wait_value(&dma_rply, D_COUNT);

  int ATHREAD_SLAVE_CORES = 64;

	// g: global, l: local, blk: block

	int num_blk_0 = (local_param.endB[0] - local_param.startB[0]) 
								* (local_param.endC[0] - local_param.startC[0]);
	int num_blk_1 = (local_param.endB[1] - local_param.startB[1]) 
								* (local_param.endC[1] - local_param.startC[1]);
	int num_blk_2 = (local_param.endB[2] - local_param.startB[2]) 
								* (local_param.endC[2] - local_param.startC[2]);
	int num_blk_3 = (local_param.endB[3] - local_param.startB[3]) 
								* (local_param.endC[3] - local_param.startC[3]);

	int g_max_num_blk = 0;
	g_max_num_blk = num_blk_0     ? num_blk_1 : num_blk_0     > num_blk_1;
	g_max_num_blk = g_max_num_blk ? num_blk_2 : g_max_num_blk > num_blk_2;
	g_max_num_blk = g_max_num_blk ? num_blk_3 : g_max_num_blk > num_blk_3;
	
	int l_max_num_blk = (g_max_num_blk + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;

	int len_a = local_param.lenA;
	int len_b = local_param.lenB;
	int len_c = local_param.lenC;

	size_t size_blk = sizeof(double) * len_a;
	int allocatable_num_blk = get_allocatable_size() / size_blk;

	int num_buf_blk = 0;
	if (l_max_num_blk < allocatable_num_blk) {
		num_buf_blk = l_max_num_blk;
	} else {
		num_buf_blk = allocatable_num_blk;
	}
	double* buffer = (double*) ldm_malloc (size_blk * num_buf_blk);

	double* arrSrc = local_param.arrSrc;
	double* arrObj = local_param.arrObj;

	int strideSrc = len_b * len_c;
	int strideObj = len_c * len_a;

	int l_num_blk, num_tile, start_b, start_c, len_blk_c;

	/* block 0 */
	l_num_blk = (num_blk_0 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[0];
	start_c = local_param.startC[0];
	len_blk_c = local_param.endC[0] - start_c;

	// arrSrc -> buffer
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if (idex_in_blk >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_0) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ia*strideSrc + ib*len_c + ic];
			}
		}
	}

	// buffer -> arrObj
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if (idex_in_blk >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_0) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				arrObj[ib*strideObj + ic*len_a + ia] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 1 */
	l_num_blk = (num_blk_1 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[1];
	start_c = local_param.startC[1];
	len_blk_c = local_param.endC[1] - start_c;

	// arrSrc -> buffer
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if (idex_in_blk >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_1) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ia*strideSrc + ib*len_c + ic];
			}
		}
	}

	// buffer -> arrObj
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if (idex_in_blk >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_1) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				arrObj[ib*strideObj + ic*len_a + ia] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 2 */
	l_num_blk = (num_blk_2 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[2];
	start_c = local_param.startC[2];
	len_blk_c = local_param.endC[2] - start_c;

	// arrSrc -> buffer
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if (idex_in_blk >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_2) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ia*strideSrc + ib*len_c + ic];
			}
		}
	}

	// buffer -> arrObj
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if (idex_in_blk >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_2) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				arrObj[ib*strideObj + ic*len_a + ia] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 3 */
	l_num_blk = (num_blk_3 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[3];
	start_c = local_param.startC[3];
	len_blk_c = local_param.endC[3] - start_c;

	// arrSrc -> buffer
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if (idex_in_blk >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_3) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ia*strideSrc + ib*len_c + ic];
			}
		}
	}

	// buffer -> arrObj
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if (idex_in_blk >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_3) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				arrObj[ib*strideObj + ic*len_a + ia] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	ldm_free (buffer, size_blk * num_buf_blk);

	return ;
}

void athread_put_halo_transpose_double_slave (struct HaloTransposeDouble3D* arg_param) {

	struct HaloTransposeDouble3D local_param;
  athread_dma_get(&local_param, arg_param, 
      sizeof(struct HaloTransposeDouble3D));

  int ATHREAD_SLAVE_CORES = 64;

	// g: global, l: local, blk: block

	int num_blk_0 = (local_param.endB[0] - local_param.startB[0]) 
								* (local_param.endC[0] - local_param.startC[0]);
	int num_blk_1 = (local_param.endB[1] - local_param.startB[1]) 
								* (local_param.endC[1] - local_param.startC[1]);
	int num_blk_2 = (local_param.endB[2] - local_param.startB[2]) 
								* (local_param.endC[2] - local_param.startC[2]);
	int num_blk_3 = (local_param.endB[3] - local_param.startB[3]) 
								* (local_param.endC[3] - local_param.startC[3]);

	int g_max_num_blk = 0;
	g_max_num_blk = num_blk_0     ? num_blk_1 : num_blk_0     > num_blk_1;
	g_max_num_blk = g_max_num_blk ? num_blk_2 : g_max_num_blk > num_blk_2;
	g_max_num_blk = g_max_num_blk ? num_blk_3 : g_max_num_blk > num_blk_3;
	
	int l_max_num_blk = (g_max_num_blk + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;

	int len_a = local_param.lenA;
	int len_b = local_param.lenB;
	int len_c = local_param.lenC;

	size_t size_blk = sizeof(double) * len_a;
	int allocatable_num_blk = get_allocatable_size() / size_blk;

	int num_buf_blk = 0;
	if (l_max_num_blk < allocatable_num_blk) {
		num_buf_blk = l_max_num_blk;
	} else {
		num_buf_blk = allocatable_num_blk;
	}
	double* buffer = (double*) ldm_malloc (size_blk * num_buf_blk);

	double* arrSrc = local_param.arrSrc;
	double* arrObj = local_param.arrObj;

	int strideSrc = len_c * len_a;
	int strideObj = len_b * len_c;

	int l_num_blk, num_tile, start_b, start_c, len_blk_c;

	/* block 0 */
	l_num_blk = (num_blk_0 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[0];
	start_c = local_param.startC[0];
	len_blk_c = local_param.endC[0] - start_c;

	// arrSrc -> buffer
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if ((idex_in_blk) >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_0) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ib*strideSrc + ic*len_a + ia];
			}
		}
	}

	// buffer -> arrObj
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if ((idex_in_blk) >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_0) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				arrObj[ia*strideObj + ib*len_c + ic] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 1 */
	l_num_blk = (num_blk_1 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[1];
	start_c = local_param.startC[1];
	len_blk_c = local_param.endC[1] - start_c;

	// arrSrc -> buffer
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if ((idex_in_blk) >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_1) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ib*strideSrc + ic*len_a + ia];
			}
		}
	}

	// buffer -> arrObj
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if ((idex_in_blk) >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_1) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				arrObj[ia*strideObj + ib*len_c + ic] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 2 */
	l_num_blk = (num_blk_2 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[2];
	start_c = local_param.startC[2];
	len_blk_c = local_param.endC[2] - start_c;

	// arrSrc -> buffer
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if ((idex_in_blk) >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_2) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ib*strideSrc + ic*len_a + ia];
			}
		}
	}

	// buffer -> arrObj
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if ((idex_in_blk) >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_2) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				arrObj[ia*strideObj + ib*len_c + ic] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	/* block 3 */
	l_num_blk = (num_blk_3 + ATHREAD_SLAVE_CORES - 1) 
			/ ATHREAD_SLAVE_CORES;
	num_tile = (l_num_blk + allocatable_num_blk - 1) / allocatable_num_blk;

	start_b = local_param.startB[3];
	start_c = local_param.startC[3];
	len_blk_c = local_param.endC[3] - start_c;

	// arrSrc -> buffer
	for (int tid = 0; tid < num_tile; ++tid) {
		for (int bid = 0; bid < allocatable_num_blk; ++bid) {
      int idex_in_blk = tid * allocatable_num_blk + bid;
			if ((idex_in_blk) >= l_num_blk) break;
			int index = athread_tid * l_num_blk + idex_in_blk;
			if (index >= num_blk_3) break;
 
			int ib = start_b + index / len_blk_c;
			if (ib >= len_b) break;
			int ic = start_c + index % len_blk_c;
			if (ic >= len_c) break;
 
			for (int ia = 0; ia < len_a; ++ia) {
				buffer[ia*num_buf_blk + idex_in_blk] = arrSrc[ib*strideSrc + ic*len_a + ia];
			}
		}
	}

	// buffer -> arrObj
	for (int ia = 0; ia < len_a; ++ia) {
		for (int tid = 0; tid < num_tile; ++tid) {
			for (int bid = 0; bid < allocatable_num_blk; ++bid) {
        int idex_in_blk = tid * allocatable_num_blk + bid;
				if ((idex_in_blk) >= l_num_blk) break;
				int index = athread_tid * l_num_blk + idex_in_blk;
				if (index >= num_blk_3) break;
				int ib = start_b + index / len_blk_c;
				if (ib >= len_b) break;
				int ic = start_c + index % len_blk_c;
				if (ic >= len_c) break;
				arrObj[ia*strideObj + ib*len_c + ic] = buffer[ia*num_buf_blk + idex_in_blk];
			}
		}
	}

	ldm_free (buffer, size_blk * num_buf_blk);

	return ;
}

#endif // __sw_slave__
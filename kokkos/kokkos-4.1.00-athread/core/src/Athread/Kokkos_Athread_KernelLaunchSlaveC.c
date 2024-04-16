#include "slave.h"

void slave_cores_reduce(void* src_addr, void* dest_addr, 
		int units, int dtype, int optype, void *buf, int buf_item) {

  athread_ssync_array();
	athread_redurt(src_addr, dest_addr, units, dtype, optype, buf, buf_item);
  // athread_ssync_array();
	return ;
}
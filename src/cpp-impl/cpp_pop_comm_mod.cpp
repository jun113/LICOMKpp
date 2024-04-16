#include "../head/def-undef.h"
#ifndef LICOM_ENABLE_FORTRAN
#include "../head/fortran_pop_comm_mod.h"

namespace CppPOPCommMod {

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &POP_communicator = pop_commmod_mp_pop_communicator_;
int &POP_myTask       = pop_commmod_mp_pop_mytask_;
int &POP_masterTask   = pop_commmod_mp_pop_mastertask_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &POP_communicator = __pop_commmod_MOD_pop_communicator;
int &POP_myTask       = __pop_commmod_MOD_pop_mytask;
int &POP_masterTask   = __pop_commmod_MOD_pop_mastertask;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppPOPCommMod

#endif // LICOM_ENABLE_FORTRAN
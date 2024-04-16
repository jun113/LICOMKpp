#include "../head/fortran_msg_mod.h"

namespace CppMsgMod {

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &nproc = msg_mod_mp_nproc_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &nproc = __msg_mod_MOD_nproc;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU

} // namespace CppMsgMod


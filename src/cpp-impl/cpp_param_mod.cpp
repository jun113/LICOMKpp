#include "../head/def-undef.h"
#include "../head/fortran_param_mod.h"

namespace CppParamMod {
#ifdef LICOM_ENABLE_FORTRAN_COMPILER_INTEL
int &ierr           = param_mod_mp_ierr_;
int &mytid          = param_mod_mp_mytid_;
                   
int &jj_start       = param_mod_mp_jj_start_;
int &jj_end         = param_mod_mp_jj_end_;
                   
int &my_task        = param_mod_mp_my_task_;
int &master_task    = param_mod_mp_master_task_;

double &am_hor      = param_mod_mp_am_hor_;
double &ah_hor      = param_mod_mp_ah_hor_;

double &am_bihar    = param_mod_mp_am_bihar_;
double &ah_bihar    = param_mod_mp_ah_bihar_;

int &num_cpl        = param_mod_mp_num_cpl_;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_INTEL

#ifdef LICOM_ENABLE_FORTRAN_COMPILER_GNU
int &ierr           = __param_mod_MOD_ierr;
int &mytid          = __param_mod_MOD_mytid;
                   
int &jj_start       = __param_mod_MOD_jj_start;
int &jj_end         = __param_mod_MOD_jj_end;
                   
int &my_task        = __param_mod_MOD_my_task;
int &master_task    = __param_mod_MOD_master_task;

double &am_hor      = __param_mod_MOD_am_hor;
double &ah_hor      = __param_mod_MOD_ah_hor;

double &am_bihar    = __param_mod_MOD_am_bihar;
double &ah_bihar    = __param_mod_MOD_ah_bihar;

int &num_cpl        = __param_mod_MOD_num_cpl;
#endif // LICOM_ENABLE_FORTRAN_COMPILER_GNU
} // namespace CppParamMod

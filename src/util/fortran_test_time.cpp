#include "../head/def-undef.h"
#include "../head/licom_test_time.hpp"

// #ifdef LICOM_ENABLE_FORTRAN 
#ifdef LICOM_ENABLE_TEST_TIME

namespace MyTest {
extern TestTime::MyTime my_time;
} // namespace MyTest

extern "C" void fortran_test_time_start_ (const char* name) {
  using MyTest::my_time;
  my_time.testTime_start (name);
  return ;
}
extern "C" void fortran_test_time_stop_ (const char* name) {
  using MyTest::my_time;
  my_time.testTime_stop (name);
  return ;
}
#endif // LICOM_ENABLE_TEST_TIME
// #endif // LICOM_ENABLE_FORTRAN
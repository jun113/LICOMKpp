#include <slave.h>
#include <math.h>

double athread_slave_sqrt (double x) {
	return sqrt(x);
}

double athread_slave_log (double x) {
	return log(x);
}
double athread_slave_atan_double (double x) {
	return atan(x);
}

double athread_slave_sin_double (double x) {
	return sin(x);
}

double athread_slave_cos_double (double x) {
	return cos(x);
}
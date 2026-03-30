#include<math.h>
#include"report.h"
#ifndef DISABLE_VALIDATION
#include"stdio.h"
#endif

void report_validation(double result, double reference, double percent_error) {
#ifndef DISABLE_VALIDATION
  double error = fabs(reference*percent_error/100);

  if(fabs(result - reference) <= error)
    printf("OK\n");
  else
    printf("NG\n");

  printf("# result = %.15e, reference = %.15e, error range = [%.15e,%.15e]\n",
         result, reference, reference - error, reference + error);
#endif
}

double get_ss_r8(double *arr, int n) {
#ifndef DISABLE_VALIDATION
  double ss = 0;
  int i;
  for(i=0; i<n; i++)
    ss += arr[i]*arr[i];
  return ss;
#else
  return 0;
#endif
}

float get_ss_r4(float *arr, int n) {
#ifndef DISABLE_VALIDATION
  double ss = 0;
  int i;
  for(i=0; i<n; i++)
    ss += arr[i]*arr[i];
  return ss;
#else
  return 0;
#endif
}

static unsigned long int next = 1;
static unsigned long int RAND_MAX_ = 32767;

static int myrand() {
  next = next * 1103515245 + 12345;
  return (unsigned int)(next/65536) % 32768;
}

double get_rand(double min, double max) {
  return ((double)myrand()/RAND_MAX_) * (max - min) + min;
}

void set_seed(unsigned int seed) {
  next = seed;
}  

void report_validation_(double *result, double *reference, double *percent_error) {
  report_validation(*result, *reference, *percent_error);
}

void get_ss_r8_(double *arr, int *n, double *result) {
  *result = get_ss_r8(arr, *n);
}

void get_ss_r4_(float *arr, int *n, float *result) {
  *result = get_ss_r4(arr, *n);
}

void get_rand_r4_(float *min, float *max, float *result) {
  *result = get_rand(*min, *max);
}

void get_rand_r8_(double *min, double *max, double *result) {
  *result = get_rand(*min, *max);
}

void set_seed_(unsigned int *seed) {
  set_seed(*seed);
}

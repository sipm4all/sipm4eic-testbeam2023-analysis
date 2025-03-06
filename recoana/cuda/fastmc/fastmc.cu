/** fastsim cuda.cu **/
#include <cuda_runtime.h>
#include <stdio.h>
#include "common.h"
//#include <iostream>

static void HandleError( cudaError_t err, const char *file, int line ) {
  if (err != cudaSuccess) {
    printf( "%s in %s at line %d\n", cudaGetErrorString( err ), file, line );
    exit( EXIT_FAILURE );
  }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

float *d_x_hits;
float *d_y_hits;
float *d_x_coords;
float *d_y_coords;
int *d_channels;

__global__ void
fastmc_process_kernel(int nhits, float *x_hits, float *y_hits, int *channels, float *x_coords, float *y_coords)
{
  int ich = blockIdx.x * blockDim.x + threadIdx.x;
  auto x = x_coords[ich];
  auto y = y_coords[ich];
  for (int ihit = 0; ihit < nhits; ++ihit) {
    if (fabsf(x_hits[ihit] - x) < 1.5 && fabsf(y_hits[ihit] - y) < 1.5)
      channels[ihit] = ich;
  }
}

void
fastmc_init(float *h_x_coords, float *h_y_coords)
{
  HANDLE_ERROR( cudaMalloc((void **)&d_x_hits, MAX_HITS * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_y_hits, MAX_HITS * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_channels, MAX_HITS * sizeof(int)) );

  HANDLE_ERROR( cudaMalloc((void **)&d_x_coords, N_CHANNELS * sizeof(float)) );
  HANDLE_ERROR( cudaMalloc((void **)&d_y_coords, N_CHANNELS * sizeof(float)) );
  
  HANDLE_ERROR( cudaMemcpy(d_x_coords, h_x_coords, N_CHANNELS * sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_y_coords, h_y_coords, N_CHANNELS * sizeof(float), cudaMemcpyHostToDevice) );
}

void
fastmc_free()
{
  cudaFree(d_x_coords);
  cudaFree(d_y_coords);
  cudaFree(d_x_hits);
  cudaFree(d_y_hits);
}

void
fastmc_process(int nhits, float_t *h_x_hits, float *h_y_hits, int *h_channels)
{
  HANDLE_ERROR( cudaMemcpy(d_x_hits, h_x_hits, nhits * sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(d_y_hits, h_y_hits, nhits * sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemset(d_channels, -1, nhits * sizeof(int)) );
  
  fastmc_process_kernel<<<N_PDUS, 256>>>(nhits, d_x_hits, d_y_hits, d_channels, d_x_coords, d_y_coords);

  HANDLE_ERROR( cudaMemcpy(h_channels, d_channels, nhits * sizeof(int), cudaMemcpyDeviceToHost) );
}

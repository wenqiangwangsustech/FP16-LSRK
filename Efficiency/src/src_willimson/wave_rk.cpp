#include "header.h"
/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/
__GLOBAL__
void wave_rk0( float * h_W, float * W, float * t_W, float * m_W, long long WStride, float DT )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		//if( i == WStride )
		//	printf( "WStride = %d\n", WStride );
		m_W[i] = W[i];
		t_W[i] = m_W[i] + beta1  * DT * h_W[i];
		W[i]   = m_W[i] + alpha2 * DT * h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk1( float * h_W, float * W, float * t_W, float * m_W, long long WStride, float DT )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		t_W[i] += beta2 * DT * h_W[i];
		W[i]    = m_W[i] + alpha3 * DT * h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk2( float * h_W, float * W, float * t_W, float * m_W, long long WStride, float DT )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		t_W[i] += beta3 * DT * h_W[i];
		W[i]    = m_W[i] + DT * h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_rk3( float * h_W, float * W, float * t_W, float * m_W, long long WStride, float DT )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		W[i] = t_W[i] +  beta4 * DT * h_W[i];
	END_CALCULATE1D( )
}

void waveRk( GRID grid, int irk, float * h_W, float * W, float * t_W, float * m_W, float DT )
{
	WAVE_RK_FUNC wave_rk[4] = { wave_rk0, wave_rk1, wave_rk2, wave_rk3 };
	long long num = grid._nx_ * grid._ny_ * grid._nz_ * WAVESIZE;

#ifdef GPU_CUDA
	dim3 threads( 1024, 1, 1 );
	dim3 blocks;
	blocks.x = ( num + threads.x - 1 ) / threads.x;
	blocks.y = 1;
	blocks.z = 1;
	wave_rk[irk]<<< blocks, threads >>>( h_W, W, t_W, m_W, num, DT );
	//wave_rk<<< blocks, threads >>>( h_W, W, t_W, m_W, num, DT );
	CHECK( cudaDeviceSynchronize( ) );
#else
	wave_rk[irk]( h_W, W, t_W, m_W, num, DT );
#endif

}




__GLOBAL__
void new_wave_rk( float * h_W, float * W, long long WStride, float B, float DT )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		W[i]   = W[i] + B * h_W[i];
	END_CALCULATE1D( )
}


void newRk( GRID grid, float * h_W, float * W, float B, float DT )
{
	long long num = grid._nx_ * grid._ny_ * grid._nz_ * WAVESIZE;

#ifdef GPU_CUDA
	dim3 threads( 1024, 1, 1 );
	dim3 blocks;
	blocks.x = ( num + threads.x - 1 ) / threads.x;
	blocks.y = 1;
	blocks.z = 1;
	new_wave_rk
	<<< blocks, threads >>>
	( h_W, W, num, B, DT );
	CHECK( cudaDeviceSynchronize( ) );
#else
	wave_rk[irk]( h_W, W, t_W, m_W, num, DT );
#endif
}

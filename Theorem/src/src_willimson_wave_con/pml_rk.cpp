/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:pml_rk.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-05
*   Discription:
*
================================================================*/

#include "header.h"
typedef void (*WAVE_RK_FUNC )( float * h_W, float * W, float * t_W, float * m_W, long long num );

__GLOBAL__
void wave_pml_rk0( float * h_W, float * W, float * t_W, float * m_W, long long WStride )
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
		t_W[i] = m_W[i] + beta1  * h_W[i];
		W[i]   = m_W[i] + alpha2 * h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_pml_rk1( float * h_W, float * W, float * t_W, float * m_W, long long WStride )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		t_W[i] += beta2 * h_W[i];
		W[i]    = m_W[i] + alpha3 * h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_pml_rk2( float * h_W, float * W, float * t_W, float * m_W, long long WStride  )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		t_W[i] += beta3 * h_W[i];
		W[i]    = m_W[i] + h_W[i];
	END_CALCULATE1D( )
}

__GLOBAL__
void wave_pml_rk3( float * h_W, float * W, float * t_W, float * m_W, long long WStride )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	CALCULATE1D( i, 0, WStride ) 
		W[i] = t_W[i] +  beta4 * h_W[i];
	END_CALCULATE1D( )
}
void pmlRk( GRID grid, MPI_BORDER border, int irk, AUX4 Aux4_1, AUX4 Aux4_2 )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	
	int nPML = grid.nPML;

	long long num = 0;
	
	WAVE_RK_FUNC pml_rk[4] = { wave_pml_rk0, wave_pml_rk1, wave_pml_rk2, wave_pml_rk3 };
	long long numx = nPML * ny * nz * WAVESIZE;
	long long numy = nPML * nx * nz * WAVESIZE;
	long long numz = nPML * nx * ny * WAVESIZE;


#ifdef GPU_CUDA
                        
	dim3 thread( 512, 1, 1 );
	dim3 blockX;

	blockX.x = ( numx + thread.x - 1 ) / thread.x;
	blockX.y = 1;
	blockX.z = 1;

	dim3 blockY;
	blockY.x = ( numy + thread.x - 1 ) / thread.x;
	blockY.y = 1;
	blockY.z = 1;

	dim3 blockZ;
	blockZ.x = ( numz + thread.x - 1 ) / thread.x;
	blockZ.y = 1;
	blockZ.z = 1;

	if ( border.isx1 )	pml_rk[irk]<<< blockX, thread >>>( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, Aux4_1.t_Aux_x.Vx, Aux4_1.m_Aux_x.Vx, numx );
	if ( border.isy1 )  pml_rk[irk]<<< blockY, thread >>>( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, Aux4_1.t_Aux_y.Vx, Aux4_1.m_Aux_y.Vx, numy );
	if ( border.isz1 )  pml_rk[irk]<<< blockZ, thread >>>( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, Aux4_1.t_Aux_z.Vx, Aux4_1.m_Aux_z.Vx, numz );
    //                                                                                                                              
	if ( border.isx2 )	pml_rk[irk]<<< blockX, thread >>>( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, Aux4_2.t_Aux_x.Vx, Aux4_2.m_Aux_x.Vx, numx );
	if ( border.isy2 )  pml_rk[irk]<<< blockY, thread >>>( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, Aux4_2.t_Aux_y.Vx, Aux4_2.m_Aux_y.Vx, numy );
	//CHECK( cudaDeviceSynchronize( ) );
                                                                                                                                        
#ifndef FREE_SURFACE
	if ( border.isz2 )  pml_rk[irk]<<< blockZ, thread >>>( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, Aux4_2.t_Aux_z.Vx, Aux4_2.m_Aux_z.Vx, numz );
#endif                                                                                                               

                                                                                                                                   
#else                                                                                                                              
                                                                                                                                   
	
                                                                                                                                   
                                                                                                                                   
	if ( border.isx1 )	pml_rk[irk]( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, Aux4_1.t_Aux_x.Vx, Aux4_1.m_Aux_x.Vx, numx );
	if ( border.isy1 )  pml_rk[irk]( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, Aux4_1.t_Aux_y.Vx, Aux4_1.m_Aux_y.Vx, numy );
	if ( border.isz1 )  pml_rk[irk]( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, Aux4_1.t_Aux_z.Vx, Aux4_1.m_Aux_z.Vx, numz );
    //                                                                                                        
	if ( border.isx2 )	pml_rk[irk]( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, Aux4_2.t_Aux_x.Vx, Aux4_2.m_Aux_x.Vx, numx );
	if ( border.isy2 )  pml_rk[irk]( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, Aux4_2.t_Aux_y.Vx, Aux4_2.m_Aux_y.Vx, numy );
#ifndef FREE_SURFACE
	if ( border.isz2 )  pml_rk[irk]( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, Aux4_2.t_Aux_z.Vx, Aux4_2.m_Aux_z.Vx, numz );
#endif                                                                                                               


#endif

}

__GLOBAL__
void new_pml_rk( float * h_W, float * W, long long WStride, float B, float DT )
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

void newPmlRk( GRID grid, MPI_BORDER border, int irk, AUX4 Aux4_1, AUX4 Aux4_2, float B, float DT )
{
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;
	
	int nPML = grid.nPML;

	long long num = 0;
	
	long long numx = nPML * ny * nz * WAVESIZE;
	long long numy = nPML * nx * nz * WAVESIZE;
	long long numz = nPML * nx * ny * WAVESIZE;


#ifdef GPU_CUDA
                        
	dim3 thread( 512, 1, 1 );
	dim3 blockX;

	blockX.x = ( numx + thread.x - 1 ) / thread.x;
	blockX.y = 1;
	blockX.z = 1;

	dim3 blockY;
	blockY.x = ( numy + thread.x - 1 ) / thread.x;
	blockY.y = 1;
	blockY.z = 1;

	dim3 blockZ;
	blockZ.x = ( numz + thread.x - 1 ) / thread.x;
	blockZ.y = 1;
	blockZ.z = 1;

	if ( border.isx1 )	new_pml_rk<<< blockX, thread >>>( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, numx, B, DT );
	if ( border.isy1 )  new_pml_rk<<< blockY, thread >>>( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, numy, B, DT );
	if ( border.isz1 )  new_pml_rk<<< blockZ, thread >>>( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, numz, B, DT );
	if ( border.isx2 )	new_pml_rk<<< blockX, thread >>>( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, numx, B, DT );
	if ( border.isy2 )  new_pml_rk<<< blockY, thread >>>( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, numy, B, DT );
                                                                                                  
#ifndef FREE_SURFACE
	if ( border.isz2 )  new_pml_rk<<< blockZ, thread >>>( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, numz, B, DT );
#endif                                                                                                               

#else                                                                                                                              
                                                                                                                                   
	if ( border.isx1 )	new_pml_rk( Aux4_1.h_Aux_x.Vx, Aux4_1.Aux_x.Vx, numx, B, DT );
	if ( border.isy1 )  new_pml_rk( Aux4_1.h_Aux_y.Vx, Aux4_1.Aux_y.Vx, numy, B, DT );
	if ( border.isz1 )  new_pml_rk( Aux4_1.h_Aux_z.Vx, Aux4_1.Aux_z.Vx, numz, B, DT );
    //                                                              
	if ( border.isx2 )	new_pml_rk( Aux4_2.h_Aux_x.Vx, Aux4_2.Aux_x.Vx, numx, B, DT );
	if ( border.isy2 )  new_pml_rk( Aux4_2.h_Aux_y.Vx, Aux4_2.Aux_y.Vx, numy, B, DT );
#ifndef FREE_SURFACE
	if ( border.isz2 )  new_pml_rk( Aux4_2.h_Aux_z.Vx, Aux4_2.Aux_z.Vx, numz, B, DT );
#endif                                                                                                               


#endif

}

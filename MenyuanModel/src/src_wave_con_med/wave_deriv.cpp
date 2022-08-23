/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:wave.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"

#ifdef PML
#define TIMES_PML_BETA_X * pml_beta_x
#define TIMES_PML_BETA_Y * pml_beta_y
#define TIMES_PML_BETA_Z * pml_beta_z
#else
#define TIMES_PML_BETA_X 
#define TIMES_PML_BETA_Y 
#define TIMES_PML_BETA_Z 
#endif


void allocWave( GRID grid, WAVE * h_W, WAVE * W, WAVE * t_W, WAVE * m_W )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	//printf( "_nx_ = %d, _ny_ = %d, _nz_ = %d\n", _nx_, _ny_, _nz_  );
	long long num = _nx_ * _ny_ * _nz_; 
		
	FLOAT * pWave = NULL;
	long long size = sizeof( FLOAT ) * num * WAVESIZE * 4;

	CHECK( Malloc( ( void ** )&pWave, size ) );
	CHECK( Memset(  pWave, 0, size ) ); 

	h_W->Vx  = pWave + 0 * num;
	h_W->Vy  = pWave + 1 * num;
	h_W->Vz  = pWave + 2 * num;
	h_W->Txx = pWave + 3 * num;
	h_W->Tyy = pWave + 4 * num;
	h_W->Tzz = pWave + 5 * num;
	h_W->Txy = pWave + 6 * num;
	h_W->Txz = pWave + 7 * num;
	h_W->Tyz = pWave + 8 * num;

	pWave 	 = pWave + 9 * num;
	

	W->Vx  = pWave + 0 * num;
	W->Vy  = pWave + 1 * num;
	W->Vz  = pWave + 2 * num;
	W->Txx = pWave + 3 * num;
	W->Tyy = pWave + 4 * num;
	W->Tzz = pWave + 5 * num;
	W->Txy = pWave + 6 * num;
	W->Txz = pWave + 7 * num;
	W->Tyz = pWave + 8 * num;

	pWave    = pWave + 9 * num;

	t_W->Vx  = pWave + 0 * num;
	t_W->Vy  = pWave + 1 * num;
	t_W->Vz  = pWave + 2 * num;
	t_W->Txx = pWave + 3 * num;
	t_W->Tyy = pWave + 4 * num;
	t_W->Tzz = pWave + 5 * num;
	t_W->Txy = pWave + 6 * num;
	t_W->Txz = pWave + 7 * num;
	t_W->Tyz = pWave + 8 * num;

	pWave 	 = pWave + 9 * num;

	m_W->Vx  = pWave + 0 * num;
	m_W->Vy  = pWave + 1 * num;
	m_W->Vz  = pWave + 2 * num;
	m_W->Txx = pWave + 3 * num;
	m_W->Tyy = pWave + 4 * num;
	m_W->Tzz = pWave + 5 * num;
	m_W->Txy = pWave + 6 * num;
	m_W->Txz = pWave + 7 * num;
	m_W->Tyz = pWave + 8 * num;



//	printf( "Alloc: h_W = %p,", h_W->Vx );
//	printf( " W = %p,",			  W->Vx );
//	printf( " t_W = %p,",		t_W->Vx );
//	printf( " m_W = %p\n",		m_W->Vx );

}


void freeWave( WAVE h_W, WAVE W, WAVE t_W, WAVE m_W )
{

	Free( h_W.Vx );
}


__GLOBAL__															
void wave_deriv( WAVE h_W, WAVE W, CONTRAVARIANT_FLOAT con, MEDIUM_FLOAT medium, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int _nx_, int _ny_, int _nz_, float rDH, int FB1, int FB2, int FB3, float DT )	
{																	

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index;

#ifdef PML
	float pml_beta_x;
	float pml_beta_y;
	float pml_beta_z;
#endif
																	
	float mu;												
	float lambda;											
	float buoyancy;											
																	
	float xi_x;   float xi_y; 	float xi_z;		
	float et_x;   float et_y; 	float et_z;		
	float zt_x;   float zt_y; 	float zt_z;		
																	
	float Txx_xi; float Tyy_xi; float Txy_xi; 	
	float Txx_et; float Tyy_et; float Txy_et; 	
	float Txx_zt; float Tyy_zt; float Txy_zt; 	
	float Txz_xi; float Tyz_xi; float Tzz_xi;	
	float Txz_et; float Tyz_et; float Tzz_et;	
	float Txz_zt; float Tyz_zt; float Tzz_zt;	
	float Vx_xi ; float Vx_et ; float Vx_zt ;	
	float Vy_xi ; float Vy_et ; float Vy_zt ;	
	float Vz_xi ; float Vz_et ; float Vz_zt ;
	float Vx1;	  float Vx2;	float Vx3;
	float Vy1;	  float Vy2;	float Vy3;
	float Vz1;	  float Vz2;	float Vz3;
	float Txx1;	  float Txx2;   float Txx3;
	float Tyy1;	  float Tyy2;	float Tyy3;
	float Tzz1;   float Tzz2;	float Tzz3;
	float Txy1;   float Txy2;	float Txy3;
	float Txz1;   float Txz2;	float Txz3;
	float Tyz1;   float Tyz2;	float Tyz3;




	float h_WVx ;  
	float h_WVy ;  
	float h_WVz ;  
	float h_WTxx;
	float h_WTyy;
	float h_WTzz;
	float h_WTxy;
	float h_WTxz;
	float h_WTyz;



	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )  
		index = INDEX( i, j, k ); 									
		mu = (float)medium.mu[index]; 											
		lambda = (float)medium.lambda[index];										
		buoyancy = (float)medium.buoyancy[index];
		buoyancy *= Crho;
		//if( index == INDEX( 100, 100, 100 ) )
		//printf( "buoyancy = %f\n", buoyancy );

#ifdef PML
		pml_beta_x = pml_beta.x[i];
		pml_beta_y = pml_beta.y[j];
		pml_beta_z = pml_beta.z[k];
#endif

		xi_x = (float)con.xi_x[index]; xi_y = (float)con.xi_y[index]; xi_z = (float)con.xi_z[index];	
		et_x = (float)con.et_x[index]; et_y = (float)con.et_y[index]; et_z = (float)con.et_z[index];	
		zt_x = (float)con.zt_x[index]; zt_y = (float)con.zt_y[index]; zt_z = (float)con.zt_z[index];	

		 Vx_xi = L( (float)W.Vx , FB1, xi ) TIMES_PML_BETA_X; 				
		 Vy_xi = L( (float)W.Vy , FB1, xi ) TIMES_PML_BETA_X; 				
		 Vz_xi = L( (float)W.Vz , FB1, xi ) TIMES_PML_BETA_X; 				
		Txx_xi = L( (float)W.Txx, FB1, xi ) TIMES_PML_BETA_X; 				
		Tyy_xi = L( (float)W.Tyy, FB1, xi ) TIMES_PML_BETA_X; 				
		Tzz_xi = L( (float)W.Tzz, FB1, xi ) TIMES_PML_BETA_X; 				
		Txy_xi = L( (float)W.Txy, FB1, xi ) TIMES_PML_BETA_X; 				
		Txz_xi = L( (float)W.Txz, FB1, xi ) TIMES_PML_BETA_X; 				
		Tyz_xi = L( (float)W.Tyz, FB1, xi ) TIMES_PML_BETA_X; 				

		 Vx_et = L( (float)W.Vx , FB2, et ) TIMES_PML_BETA_Y;				
		 Vy_et = L( (float)W.Vy , FB2, et ) TIMES_PML_BETA_Y;				
		 Vz_et = L( (float)W.Vz , FB2, et ) TIMES_PML_BETA_Y;				
		Txx_et = L( (float)W.Txx, FB2, et ) TIMES_PML_BETA_Y;				
		Tyy_et = L( (float)W.Tyy, FB2, et ) TIMES_PML_BETA_Y;				
		Tzz_et = L( (float)W.Tzz, FB2, et ) TIMES_PML_BETA_Y;				
		Txy_et = L( (float)W.Txy, FB2, et ) TIMES_PML_BETA_Y;				
		Txz_et = L( (float)W.Txz, FB2, et ) TIMES_PML_BETA_Y;				
		Tyz_et = L( (float)W.Tyz, FB2, et ) TIMES_PML_BETA_Y;				

   		 Vx_zt = L( (float)W.Vx , FB3, zt ) TIMES_PML_BETA_Z;				
   		 Vy_zt = L( (float)W.Vy , FB3, zt ) TIMES_PML_BETA_Z;				
   		 Vz_zt = L( (float)W.Vz , FB3, zt ) TIMES_PML_BETA_Z;				
  		Txx_zt = L( (float)W.Txx, FB3, zt ) TIMES_PML_BETA_Z;				
  		Tyy_zt = L( (float)W.Tyy, FB3, zt ) TIMES_PML_BETA_Z;				
  		Tzz_zt = L( (float)W.Tzz, FB3, zt ) TIMES_PML_BETA_Z;				
  		Txy_zt = L( (float)W.Txy, FB3, zt ) TIMES_PML_BETA_Z;				
  		Txz_zt = L( (float)W.Txz, FB3, zt ) TIMES_PML_BETA_Z;				
  		Tyz_zt = L( (float)W.Tyz, FB3, zt ) TIMES_PML_BETA_Z;

		Vx1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txx_xi, Txy_xi, Txz_xi ) * buoyancy;
		Vx2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txx_et, Txy_et, Txz_et ) * buoyancy;
		Vx3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txx_zt, Txy_zt, Txz_zt ) * buoyancy;
		Vy1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txy_xi, Tyy_xi, Tyz_xi ) * buoyancy;
		Vy2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txy_et, Tyy_et, Tyz_et ) * buoyancy;
		Vy3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txy_zt, Tyy_zt, Tyz_zt ) * buoyancy;
		Vz1  = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Txz_xi, Tyz_xi, Tzz_xi ) * buoyancy;
		Vz2  = DOT_PRODUCT3D( et_x, et_y, et_z, Txz_et, Tyz_et, Tzz_et ) * buoyancy;
		Vz3  = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Txz_zt, Tyz_zt, Tzz_zt ) * buoyancy;

		Txx1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_x * Vx_xi );
		Txx2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_x * Vx_et );
		Txx3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_x * Vx_zt );
		Tyy1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_y * Vy_xi );
		Tyy2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_y * Vy_et );
		Tyy3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_y * Vy_zt );
		Tzz1 = DOT_PRODUCT3D( xi_x, xi_y, xi_z, Vx_xi, Vy_xi, Vz_xi ) * lambda + 2.0f * mu * ( xi_z * Vz_xi );
		Tzz2 = DOT_PRODUCT3D( et_x, et_y, et_z, Vx_et, Vy_et, Vz_et ) * lambda + 2.0f * mu * ( et_z * Vz_et );
		Tzz3 = DOT_PRODUCT3D( zt_x, zt_y, zt_z, Vx_zt, Vy_zt, Vz_zt ) * lambda + 2.0f * mu * ( zt_z * Vz_zt );

		Txy1 = DOT_PRODUCT2D( xi_y, xi_x, Vx_xi, Vy_xi ) * mu;
		Txy2 = DOT_PRODUCT2D( et_y, et_x, Vx_et, Vy_et ) * mu;
		Txy3 = DOT_PRODUCT2D( zt_y, zt_x, Vx_zt, Vy_zt ) * mu;
		Txz1 = DOT_PRODUCT2D( xi_z, xi_x, Vx_xi, Vz_xi ) * mu;
		Txz2 = DOT_PRODUCT2D( et_z, et_x, Vx_et, Vz_et ) * mu;
		Txz3 = DOT_PRODUCT2D( zt_z, zt_x, Vx_zt, Vz_zt ) * mu;
		Tyz1 = DOT_PRODUCT2D( xi_z, xi_y, Vy_xi, Vz_xi ) * mu;
		Tyz2 = DOT_PRODUCT2D( et_z, et_y, Vy_et, Vz_et ) * mu;
		Tyz3 = DOT_PRODUCT2D( zt_z, zt_y, Vy_zt, Vz_zt ) * mu;

		h_WVx  	= ( Vx1  + Vx2  + Vx3  ) * DT;							
		h_WVy  	= ( Vy1  + Vy2  + Vy3  ) * DT;							
		h_WVz  	= ( Vz1  + Vz2  + Vz3  ) * DT;							
		h_WTxx 	= ( Txx1 + Txx2 + Txx3 ) * DT;						
		h_WTyy 	= ( Tyy1 + Tyy2 + Tyy3 ) * DT;						
		h_WTzz 	= ( Tzz1 + Tzz2 + Tzz3 ) * DT;						
		h_WTxy 	= ( Txy1 + Txy2 + Txy3 ) * DT;						
		h_WTxz 	= ( Txz1 + Txz2 + Txz3 ) * DT;						
		h_WTyz 	= ( Tyz1 + Tyz2 + Tyz3 ) * DT;						

		h_W.Vx [index] 	= h_WVx ;							
		h_W.Vy [index] 	= h_WVy ;							
		h_W.Vz [index] 	= h_WVz ;							
		h_W.Txx[index] 	= h_WTxx;						
		h_W.Tyy[index] 	= h_WTyy;						
		h_W.Tzz[index] 	= h_WTzz;						
		h_W.Txy[index] 	= h_WTxy;						
		h_W.Txz[index] 	= h_WTxz;						
		h_W.Tyz[index] 	= h_WTyz;						
																				
																	
	END_CALCULATE3D( )												
}			


void waveDeriv( GRID grid, WAVE h_W, WAVE W, CONTRAVARIANT_FLOAT con, MEDIUM_FLOAT medium, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3, float DT )	
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	float rDH = grid.rDH;

	//printf( "Deriv: h_W = %p,", h_W.Vx );
	//printf( " W = %p\n",	  W.Vx );
	//printf( " t_W = %p,",	t_W );
	//printf( " m_W = %p\n",	m_W );

#ifdef GPU_CUDA

	int nx = _nx_ - 2 * HALO;
	int ny = _ny_ - 2 * HALO;
	int nz = _nz_ - 2 * HALO;

	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;

	wave_deriv <<< blocks, threads >>>
	( h_W, W, con, medium, 
#ifdef PML
	pml_beta,
#endif  //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );	

	CHECK( cudaDeviceSynchronize( ) );

#else   //GPU_CUDA

	wave_deriv 
	( h_W, W, con, medium, 
#ifdef PML
	pml_beta,
#endif //PML
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );	


#endif //GPU_CUDA
}





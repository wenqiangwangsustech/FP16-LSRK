/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:propagate.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-03
*   Discription:
*
================================================================*/

#include "header.h"
#ifdef PML
#define times_pml_beta_x * pml_beta_x
#define times_pml_beta_y * pml_beta_y
#define times_pml_beta_z * pml_beta_z
#else
#define times_pml_beta_x 
#define times_pml_beta_y 
#define times_pml_beta_z 
#endif

__GLOBAL__	
void free_surface_deriv(
WAVE h_W, WAVE W, CONTRAVARIANT_FLOAT con, 	MEDIUM_FLOAT medium, FLOAT * Jac, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY, 
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
	int k = threadIdx.z + blockIdx.z * blockDim.z + _nz - HALO;			
#else
	int i = HALO;
	int j = HALO;
	int k = _nz - HALO;
#endif

	long long index;


	float mu = 0.0f;												
	float lambda = 0.0f;											
	float buoyancy = 0.0f;


#ifdef PML
	float pml_beta_x = 0.0f;
	float pml_beta_y = 0.0f;
#endif
	
																	
	float xi_x = 0.0f; 	float xi_y = 0.0f; 	float xi_z = 0.0f;		
	float et_x = 0.0f; 	float et_y = 0.0f; 	float et_z = 0.0f;		
	float zt_x = 0.0f; 	float zt_y = 0.0f; 	float zt_z = 0.0f;		
																	
	float Vx_xi  = 0.0f; float Vx_et  = 0.0f; float Vx_zt  = 0.0f;	
	float Vy_xi  = 0.0f; float Vy_et  = 0.0f; float Vy_zt  = 0.0f;	
	float Vz_xi  = 0.0f; float Vz_et  = 0.0f; float Vz_zt  = 0.0f;	
																																		
	float Jinv = 0.0f;																	
	float jacb = 0.0f;
	float J_T1x[7] = { 0.0f }; 	float J_T2x[7] = { 0.0f }; 	float J_T3x[7] = { 0.0f };	
	float J_T1y[7] = { 0.0f }; 	float J_T2y[7] = { 0.0f }; 	float J_T3y[7] = { 0.0f };	
	float J_T1z[7] = { 0.0f }; 	float J_T2z[7] = { 0.0f }; 	float J_T3z[7] = { 0.0f };	


	float Txx1 = 0.0; float Txx2 = 0.0f; float Txx3 = 0.0f;
	float Tyy1 = 0.0; float Tyy2 = 0.0f; float Tyy3 = 0.0f;
	float Tzz1 = 0.0; float Tzz2 = 0.0f; float Tzz3 = 0.0f;
	float Txy1 = 0.0; float Txy2 = 0.0f; float Txy3 = 0.0f;
	float Txz1 = 0.0; float Txz2 = 0.0f; float Txz3 = 0.0f;
	float Tyz1 = 0.0; float Tyz2 = 0.0f; float Tyz3 = 0.0f;	



	float h_WVx ;  
	float h_WVy ;  
	float h_WVz ;  
	float h_WTxx;
	float h_WTyy;
	float h_WTzz;
	float h_WTxy;
	float h_WTxz;
	float h_WTyz;


																
	int l = 0;													
	int k_s = 0;/*relative index on the surface*/				
	int pos = 0;												
	
	int indexOnSurf;


	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, _nz - HALO, _nz )		
		index = INDEX( i, j, k ); 								
		mu = (float)medium.mu[index]; 									
		lambda = (float)medium.lambda[index];							
		buoyancy = (float)medium.buoyancy[index];
		buoyancy *= Crho;


#ifdef PML
		pml_beta_x = pml_beta.x[i];
		pml_beta_y = pml_beta.y[j];
#endif
																
		k_s = ( _nz - 1 ) - k + HALO;/*relative index on the surface*/
		for ( l = 0; l <= ( 2 * HALO ); l ++ )								
		{														
			pos = INDEX( i + ( l - HALO ), j, k );					
			xi_x = (float)con.xi_x[pos]; 	xi_y = (float)con.xi_y[pos]; 	xi_z = (float)con.xi_z[pos];
			jacb = Jac[pos];
			J_T1x[l] = ( xi_x * (float)W.Txx[pos] + xi_y * (float)W.Txy[pos] + xi_z * (float)W.Txz[pos] ) * jacb;
			J_T2x[l] = ( xi_x * (float)W.Txy[pos] + xi_y * (float)W.Tyy[pos] + xi_z * (float)W.Tyz[pos] ) * jacb;
			J_T3x[l] = ( xi_x * (float)W.Txz[pos] + xi_y * (float)W.Tyz[pos] + xi_z * (float)W.Tzz[pos] ) * jacb;
												
			pos = INDEX( i, j + ( l - HALO ), k );	
			et_x = (float)con.et_x[pos]; 	et_y = (float)con.et_y[pos]; 	et_z = (float)con.et_z[pos];
			jacb = Jac[pos];
			J_T1y[l] = ( et_x * (float)W.Txx[pos] + et_y * (float)W.Txy[pos] + et_z * (float)W.Txz[pos] ) * jacb;
			J_T2y[l] = ( et_x * (float)W.Txy[pos] + et_y * (float)W.Tyy[pos] + et_z * (float)W.Tyz[pos] ) * jacb;
			J_T3y[l] = ( et_x * (float)W.Txz[pos] + et_y * (float)W.Tyz[pos] + et_z * (float)W.Tzz[pos] ) * jacb;
		}										
		for ( l = 0; l < k_s; l ++ )			
		{										
			pos = INDEX( i, j, k + ( l - HALO ) );
			zt_x = (float)con.zt_x[pos]; 	zt_y = (float)con.zt_y[pos]; 	zt_z = (float)con.zt_z[pos];
			jacb = Jac[pos];
			J_T1z[l] = ( zt_x * (float)W.Txx[pos] + zt_y * (float)W.Txy[pos] + zt_z * (float)W.Txz[pos] ) * jacb;
			J_T2z[l] = ( zt_x * (float)W.Txy[pos] + zt_y * (float)W.Tyy[pos] + zt_z * (float)W.Tyz[pos] ) * jacb;
			J_T3z[l] = ( zt_x * (float)W.Txz[pos] + zt_y * (float)W.Tyz[pos] + zt_z * (float)W.Tzz[pos] ) * jacb;
		}									
		/*The T on the surface is 0.*/		
		J_T1z[k_s] = 0.0f;					
		J_T2z[k_s] = 0.0f;					
		J_T3z[k_s] = 0.0f;					
		for ( l = k_s + 1; l <= 2 * HALO; l ++ )	
		{									
			J_T1z[l] = - J_T1z[2 * k_s - l];
			J_T2z[l] = - J_T2z[2 * k_s - l];
			J_T3z[l] = - J_T3z[2 * k_s - l];
		}									
		jacb = Jac[index];
		Jinv = 1.0f / jacb;
		
		h_WVx = buoyancy * Jinv * ( L_J_T( J_T1x, FB1 ) times_pml_beta_x
		        				  + L_J_T( J_T1y, FB2 ) times_pml_beta_y
		        				  + L_J_T( J_T1z, FB3 ) ); 			
		h_WVy = buoyancy * Jinv * ( L_J_T( J_T2x, FB1 ) times_pml_beta_x
		        				  + L_J_T( J_T2y, FB2 ) times_pml_beta_y
		        				  + L_J_T( J_T2z, FB3 ) ); 			
		h_WVz = buoyancy * Jinv * ( L_J_T( J_T3x, FB1 ) times_pml_beta_x
										  + L_J_T( J_T3y, FB2 ) times_pml_beta_y
										  + L_J_T( J_T3z, FB3 ) ); 			

		Vx_xi = L( (float)W.Vx, FB1, xi ) times_pml_beta_x;	 Vx_et = L( (float)W.Vx, FB2, et ) times_pml_beta_y;
		Vy_xi = L( (float)W.Vy, FB1, xi ) times_pml_beta_x;	 Vy_et = L( (float)W.Vy, FB2, et ) times_pml_beta_y;
		Vz_xi = L( (float)W.Vz, FB1, xi ) times_pml_beta_x;	 Vz_et = L( (float)W.Vz, FB2, et ) times_pml_beta_y;

		xi_x = (float)con.xi_x[index]; 	xi_y = (float)con.xi_y[index]; 	xi_z = (float)con.xi_z[index];
		et_x = (float)con.et_x[index]; 	et_y = (float)con.et_y[index]; 	et_z = (float)con.et_z[index];
		zt_x = (float)con.zt_x[index]; 	zt_y = (float)con.zt_y[index]; 	zt_z = (float)con.zt_z[index];

//=======================================================
//When change the HALO, BE CAREFUL!!!!
//=======================================================
											
		if ( k == _nz - 1 )					
		{ 									
			indexOnSurf = INDEX( i, j, 0 );	
			Vx_zt = DOT_PRODUCT3D( _rDZ_DX.M11[indexOnSurf], _rDZ_DX.M12[indexOnSurf], _rDZ_DX.M13[indexOnSurf], Vx_xi, Vy_xi, Vz_xi ) 	
				  + DOT_PRODUCT3D( _rDZ_DY.M11[indexOnSurf], _rDZ_DY.M12[indexOnSurf], _rDZ_DY.M13[indexOnSurf], Vx_et, Vy_et, Vz_et );	
			Vy_zt = DOT_PRODUCT3D( _rDZ_DX.M21[indexOnSurf], _rDZ_DX.M22[indexOnSurf], _rDZ_DX.M23[indexOnSurf], Vx_xi, Vy_xi, Vz_xi ) 	
				  + DOT_PRODUCT3D( _rDZ_DY.M21[indexOnSurf], _rDZ_DY.M22[indexOnSurf], _rDZ_DY.M23[indexOnSurf], Vx_et, Vy_et, Vz_et );	
			Vz_zt = DOT_PRODUCT3D( _rDZ_DX.M31[indexOnSurf], _rDZ_DX.M32[indexOnSurf], _rDZ_DX.M33[indexOnSurf], Vx_xi, Vy_xi, Vz_xi ) 	
				  + DOT_PRODUCT3D( _rDZ_DY.M31[indexOnSurf], _rDZ_DY.M32[indexOnSurf], _rDZ_DY.M33[indexOnSurf], Vx_et, Vy_et, Vz_et );	
		} 								
		if ( k == _nz - 2 ) 				
		{								
			Vx_zt =	L2( (float)W.Vx, FB3, zt );
			Vy_zt =	L2( (float)W.Vy, FB3, zt );
			Vz_zt =	L2( (float)W.Vz, FB3, zt );
		}								
		if ( k == _nz - 3 ) 				
		{								
			Vx_zt =	L3( (float)W.Vx, FB3, zt );
			Vy_zt =	L3( (float)W.Vy, FB3, zt );
			Vz_zt =	L3( (float)W.Vz, FB3, zt );
		}								
		
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
											
		h_WTxx = Txx1 + Txx2 + Txx3;
		h_WTyy = Tyy1 + Tyy2 + Tyy3;
		h_WTzz = Tzz1 + Tzz2 + Tzz3;
		h_WTxy = Txy1 + Txy2 + Txy3;
		h_WTxz = Txz1 + Txz2 + Txz3;
		h_WTyz = Tyz1 + Tyz2 + Tyz3;

		h_W.Vx [index] 	= h_WVx  * DT;							
		h_W.Vy [index] 	= h_WVy  * DT;							
		h_W.Vz [index] 	= h_WVz  * DT;							
		h_W.Txx[index] 	= h_WTxx * DT;						
		h_W.Tyy[index] 	= h_WTyy * DT;						
		h_W.Tzz[index] 	= h_WTzz * DT;						
		h_W.Txy[index] 	= h_WTxy * DT;						
		h_W.Txz[index] 	= h_WTxz * DT;						
		h_W.Tyz[index] 	= h_WTyz * DT;						
	END_CALCULATE3D( )																													
}			


void freeSurfaceDeriv( 
	GRID grid, WAVE h_W, WAVE W, CONTRAVARIANT_FLOAT con, 	
	MEDIUM_FLOAT medium, FLOAT * Jac, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3, float DT )							
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	float rDH = grid.rDH;

#ifdef GPU_CUDA
	int nx = _nx_ - 2 * HALO;
	int ny = _ny_ - 2 * HALO;
	//int nz = _nz_ - 2 * HALO;

	dim3 threads( 32, 4, 1);
	dim3 blocks;
	
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = HALO / threads.z;

	free_surface_deriv <<< blocks, threads >>>
	( h_W, W, con, medium, Jac, _rDZ_DX, _rDZ_DY, 
#ifdef PML
	pml_beta,
#endif
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );

#else

	free_surface_deriv
	( h_W, W, con, medium, Jac, _rDZ_DX, _rDZ_DY, 
#ifdef PML
	pml_beta,
#endif
	_nx_, _ny_, _nz_, rDH, FB1, FB2, FB3, DT );


#endif //GPU_CUDA
}


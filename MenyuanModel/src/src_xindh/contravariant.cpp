/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:matrix.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-02
*   Discription:
*
================================================================*/
#include "header.h"



void allocContravariant( GRID grid, CONTRAVARIANT * con )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

		
	float * pContravariant = NULL;
	long long size = sizeof( float ) * num * CONTRASIZE;

	CHECK( Malloc( ( void ** )&pContravariant, size ) );
	CHECK( Memset(  pContravariant, 0, size ) ); 

	con->xi_x = pContravariant + num * 0;
	con->xi_y = pContravariant + num * 1;
	con->xi_z = pContravariant + num * 2;
                                                  
	con->et_x = pContravariant + num * 3;
	con->et_y = pContravariant + num * 4;
	con->et_z = pContravariant + num * 5;
                                                  
	con->zt_x = pContravariant + num * 6;
	con->zt_y = pContravariant + num * 7;
	con->zt_z = pContravariant + num * 8;

}

void freeContravariant( CONTRAVARIANT con )
{	
	Free( con.xi_x );
	//con.xi_x = NULL;
}

void allocContravariant_FLOAT( GRID grid, CONTRAVARIANT_FLOAT * con )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

		
	FLOAT * pContravariant = NULL;
	long long size = sizeof( FLOAT ) * num * CONTRASIZE;

	CHECK( Malloc( ( void ** )&pContravariant, size ) );
	CHECK( Memset(  pContravariant, 0, size ) ); 

	con->xi_x = pContravariant + num * 0;
	con->xi_y = pContravariant + num * 1;
	con->xi_z = pContravariant + num * 2;
                                                  
	con->et_x = pContravariant + num * 3;
	con->et_y = pContravariant + num * 4;
	con->et_z = pContravariant + num * 5;
                                                  
	con->zt_x = pContravariant + num * 6;
	con->zt_y = pContravariant + num * 7;
	con->zt_z = pContravariant + num * 8;

}

void freeContravariant_FLOAT( CONTRAVARIANT_FLOAT con )
{	
	Free( con.xi_x );
	//con.xi_x = NULL;
}

void allocJac( GRID grid, float ** Jac )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

	float * pJac = NULL;
	long long size = sizeof( float ) * num;

	Malloc( ( void ** )&pJac, size );
	Memset(  pJac, 0, size ); 
	
	*Jac = pJac;

}


void freeJac( float * Jac )
{
	Free( Jac );
}

void allocJac_FLOAT( GRID grid, FLOAT ** Jac )
{

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;
	
	long long num = _nx_ * _ny_ * _nz_; 

	FLOAT * pJac = NULL;
	long long size = sizeof( FLOAT ) * num;

	Malloc( ( void ** )&pJac, size );
	Memset(  pJac, 0, size ); 
	
	*Jac = pJac;

}


void freeJac_FLOAT( FLOAT * Jac )
{
	Free( Jac );
}

//When change the fast axises:
/*
 * =============================================
 *             BE careful!!!!!!!!!!!!
 * =============================================
*/

__GLOBAL__
void solve_con_jac( CONTRAVARIANT con, COORD coord, float * Jac, int _nx_, int _ny_, int _nz_, float rDH )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + HALO;
#else
	int i = HALO;
	int j = HALO;
	int k = HALO;
#endif
	long long index = 0;

	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;


	//float rDH = 1.0f / DH;

	float JacobiInv = 0.0f;

	float x_xi = 0.0f, x_et = 0.0f, x_zt = 0.0f;
	float y_xi = 0.0f, y_et = 0.0f, y_zt = 0.0f;
	float z_xi = 0.0f, z_et = 0.0f, z_zt = 0.0f;

	float xi_x = 0.0f, et_x = 0.0f, zt_x = 0.0f;
	float xi_y = 0.0f, et_y = 0.0f, zt_y = 0.0f;
	float xi_z = 0.0f, et_z = 0.0f, zt_z = 0.0f;

	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, HALO, _nz )
		index = INDEX( i, j, k );
		x_xi = 0.5f * ( L( coord.x, 1,  xi ) + L( coord.x, -1, xi ) );	x_et = 0.5f * ( L( coord.x, 1, et ) + L( coord.x, -1, et ) );	x_zt = 0.5f * ( L( coord.x, 1, zt ) + L( coord.x, -1, zt ) );
		y_xi = 0.5f * ( L( coord.y, 1,  xi ) + L( coord.y, -1, xi ) );	y_et = 0.5f * ( L( coord.y, 1, et ) + L( coord.y, -1, et ) );	y_zt = 0.5f * ( L( coord.y, 1, zt ) + L( coord.y, -1, zt ) );
		z_xi = 0.5f * ( L( coord.z, 1,  xi ) + L( coord.z, -1, xi ) );	z_et = 0.5f * ( L( coord.z, 1, et ) + L( coord.z, -1, et ) );	z_zt = 0.5f * ( L( coord.z, 1, zt ) + L( coord.z, -1, zt ) );

		Jac[index] = x_xi * y_et * z_zt + x_et * y_zt * z_xi + x_zt * y_xi * z_et - x_zt * y_et * z_xi - x_et * y_xi * z_zt - x_xi * z_et * y_zt;
	  	JacobiInv = 1.0f / Jac[index];

	  	xi_x = ( y_et * z_zt - y_zt * z_et ) * JacobiInv;	xi_y = ( x_zt * z_et - x_et * z_zt ) * JacobiInv;	xi_z = ( x_et * y_zt - x_zt * y_et ) * JacobiInv;
	  	et_x = ( y_zt * z_xi - y_xi * z_zt ) * JacobiInv;	et_y = ( x_xi * z_zt - x_zt * z_xi ) * JacobiInv;	et_z = ( x_zt * y_xi - x_xi * y_zt ) * JacobiInv;
	  	zt_x = ( y_xi * z_et - y_et * z_xi ) * JacobiInv;	zt_y = ( x_et * z_xi - x_xi * z_et ) * JacobiInv;	zt_z = ( x_xi * y_et - x_et * y_xi ) * JacobiInv;

		con.xi_x[index] = xi_x; con.xi_y[index] = xi_y; con.xi_z[index] = xi_z;
		con.et_x[index] = et_x; con.et_y[index] = et_y; con.et_z[index] = et_z;
		con.zt_x[index] = zt_x; con.zt_y[index] = zt_y; con.zt_z[index] = zt_z;

	END_CALCULATE3D( )
}

//When change the fast axises:
/*
 * =============================================
 *             BE careful!!!!!!!!!!!!
 * =============================================
*/
void allocMat3x3( GRID grid, Mat3x3 * _rDZ_DX, Mat3x3 * _rDZ_DY )
{
	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	
	long long num = _nx_ * _ny_; 
		
	float * pSurf = NULL;
	long long size = sizeof( float ) * num * 2 * 9;

	Malloc( ( void ** )&pSurf, size );
	Memset(  pSurf, 0, size ); 
	
	_rDZ_DX->M11 = pSurf + 0 * num; 
	_rDZ_DX->M12 = pSurf + 1 * num; 
	_rDZ_DX->M13 = pSurf + 2 * num;
	_rDZ_DX->M21 = pSurf + 3 * num; 
	_rDZ_DX->M22 = pSurf + 4 * num; 
	_rDZ_DX->M23 = pSurf + 5 * num;
	_rDZ_DX->M31 = pSurf + 6 * num; 
	_rDZ_DX->M32 = pSurf + 7 * num; 
	_rDZ_DX->M33 = pSurf + 8 * num;

	pSurf		 = pSurf + 9 * num; 

	_rDZ_DY->M11 = pSurf + 0 * num; 
	_rDZ_DY->M12 = pSurf + 1 * num; 
	_rDZ_DY->M13 = pSurf + 2 * num;
	_rDZ_DY->M21 = pSurf + 3 * num; 
	_rDZ_DY->M22 = pSurf + 4 * num; 
	_rDZ_DY->M23 = pSurf + 5 * num;
	_rDZ_DY->M31 = pSurf + 6 * num; 
	_rDZ_DY->M32 = pSurf + 7 * num; 
	_rDZ_DY->M33 = pSurf + 8 * num;

}

void freeMat3x3( Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY )
{
	Free( _rDZ_DX.M11 );
}


__GLOBAL__
void solve_coordinate_on_free_surface( 
	CONTRAVARIANT con, COORD coord, float * Jac,
	MEDIUM medium, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY, 
	int _nx_, int _ny_, int _nz_ ) 
{
	int _nx = _nx_ - HALO;
	int _ny = _ny_ - HALO;
	int _nz = _nz_ - HALO;
	

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x + HALO; 
	int j = threadIdx.y + blockIdx.y * blockDim.y + HALO;
	int k = threadIdx.z + blockIdx.z * blockDim.z + _nz - 1;
#else
	int i = HALO;
	int j = HALO;
	int k = _nz - 1;
#endif

	long long index = 0;
	long long pos = 0;

	float DZ[9];
	float DX[9];
	float DY[9];
	float rDZ[9]; //the inverse matrix of DZ
	float DZ_det;
	float chi = 0.0f;
	float lambda = 0.0f;
	float mu = 0.0f;
	float lam_mu = 0.0f;
	


	int indexOnSurf = 0;

	float xi_x = 0.0f, et_x = 0.0f, zt_x = 0.0f;
	float xi_y = 0.0f, et_y = 0.0f, zt_y = 0.0f;
	float xi_z = 0.0f, et_z = 0.0f, zt_z = 0.0f;	


	CALCULATE3D( i, j, k, HALO, _nx, HALO, _ny, _nz - 1, _nz_ )
		index = INDEX( i, j, k );

		if ( k == _nz - 1 )
		{	
			 xi_x = con.xi_x[index];  xi_y = con.xi_y[index];  xi_z = con.xi_z[index];
			 et_x = con.et_x[index];  et_y = con.et_y[index];  et_z = con.et_z[index];
			 zt_x = con.zt_x[index];  zt_y = con.zt_y[index];  zt_z = con.zt_z[index];

			indexOnSurf = INDEX( i, j, 0 );
			lambda = medium.lambda[index];
			mu = medium.mu[index];
			chi = 2.0f * lambda + mu;//希腊字母
			lam_mu = lambda + mu;
			/****************
			---				  ---
			| DZ[0] DZ[1] DZ[2] |
			| DZ[3] DZ[4] DZ[5] |
			| DZ[6] DZ[7] DZ[8] |
			---				  ---
			*****************/
			DZ[0] = chi * zt_x * zt_x + mu * ( zt_y * zt_y + zt_z * zt_z );	DZ[1] = lam_mu * zt_x * zt_y; 									DZ[2] = lam_mu * zt_x * zt_z;
			DZ[3] = DZ[1];													DZ[4] = chi * zt_y * zt_y + mu * ( zt_x * zt_x + zt_z * zt_z );	DZ[5] = lam_mu * zt_y * zt_z;
			DZ[6] = DZ[2];													DZ[7] = DZ[5];													DZ[8] = chi * zt_z * zt_z + mu * ( zt_x * zt_x + zt_y * zt_y );

			DZ_det 	= DZ[0] * DZ[4] * DZ[8] 
					+ DZ[1] * DZ[5] * DZ[6] 
					+ DZ[2] * DZ[7] * DZ[3] 
					- DZ[2] * DZ[4] * DZ[6] 
					- DZ[1] * DZ[3] * DZ[8] 
					- DZ[0] * DZ[7] * DZ[5];

			DX[0] = chi * zt_x * xi_x + mu * ( zt_y * xi_y + zt_z * xi_z );		DX[1] = lambda * zt_x * xi_y + mu * zt_y * xi_x;					DX[2] = lambda * zt_x * xi_z + mu * zt_z * xi_x;
			DX[3] = lambda * zt_y * xi_x + mu * zt_x * xi_y;					DX[4] = chi * zt_y * xi_y + mu * ( zt_x * xi_x + zt_z * xi_z );		DX[5] = lambda * zt_y * xi_z + mu * zt_z * xi_y;
			DX[6] = lambda * zt_z * xi_x + mu * zt_x * xi_z;					DX[7] = lambda * zt_z * xi_y + mu * zt_y * xi_z;					DX[8] = chi * zt_z * xi_z + mu * ( zt_x * xi_x + zt_y * xi_y );

			DY[0] = chi * zt_x * et_x + mu * ( zt_y * et_y + zt_z * et_z );		DY[1] = lambda * zt_x * et_y + mu * zt_y * et_x;					DY[2] = lambda * zt_x * et_z + mu * zt_z * et_x;
			DY[3] = lambda * zt_y * et_x + mu * zt_x * et_y;					DY[4] = chi * zt_y * et_y + mu * ( zt_x * et_x + zt_z * et_z );		DY[5] = lambda * zt_y * et_z + mu * zt_z * et_y;
			DY[6] = lambda * zt_z * et_x + mu * zt_x * et_z;					DY[7] = lambda * zt_z * et_y + mu * zt_y * et_z;					DY[8] = chi * zt_z * et_z + mu * ( zt_x * et_x + zt_y * et_y );

			rDZ[0] = (   DZ[4] * DZ[8] - DZ[5] * DZ[7] ) / DZ_det; 
			rDZ[1] = ( - DZ[3] * DZ[8] + DZ[5] * DZ[6] ) / DZ_det; 
			rDZ[2] = (   DZ[3] * DZ[7] - DZ[4] * DZ[6] ) / DZ_det; 
			rDZ[3] = ( - DZ[1] * DZ[8] + DZ[2] * DZ[7] ) / DZ_det; 
			rDZ[4] = (   DZ[0] * DZ[8] - DZ[2] * DZ[6] ) / DZ_det; 
			rDZ[5] = ( - DZ[0] * DZ[7] + DZ[1] * DZ[6] ) / DZ_det; 
			rDZ[6] = (   DZ[1] * DZ[5] - DZ[2] * DZ[4] ) / DZ_det; 
			rDZ[7] = ( - DZ[0] * DZ[5] + DZ[2] * DZ[3] ) / DZ_det; 
			rDZ[8] = (   DZ[0] * DZ[4] - DZ[1] * DZ[3] ) / DZ_det;
			
			_rDZ_DX.M11[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[0], DX[3], DX[6] );
			_rDZ_DX.M12[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[1], DX[4], DX[7] );
			_rDZ_DX.M13[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DX[2], DX[5], DX[8] );

			_rDZ_DX.M21[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[0], DX[3], DX[6] );
			_rDZ_DX.M22[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[1], DX[4], DX[7] );
			_rDZ_DX.M23[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DX[2], DX[5], DX[8] );

			_rDZ_DX.M31[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[0], DX[3], DX[6] );
			_rDZ_DX.M32[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[1], DX[4], DX[7] );
			_rDZ_DX.M33[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DX[2], DX[5], DX[8] );

			_rDZ_DY.M11[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[0], DY[3], DY[6] );
			_rDZ_DY.M12[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[1], DY[4], DY[7] );
			_rDZ_DY.M13[indexOnSurf] = - DOT_PRODUCT3D( rDZ[0], rDZ[1], rDZ[2], DY[2], DY[5], DY[8] );

			_rDZ_DY.M21[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[0], DY[3], DY[6] );
			_rDZ_DY.M22[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[1], DY[4], DY[7] );
			_rDZ_DY.M23[indexOnSurf] = - DOT_PRODUCT3D( rDZ[3], rDZ[4], rDZ[5], DY[2], DY[5], DY[8] );
			
			_rDZ_DY.M31[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[0], DY[3], DY[6] );
			_rDZ_DY.M32[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[1], DY[4], DY[7] );
			_rDZ_DY.M33[indexOnSurf] = - DOT_PRODUCT3D( rDZ[6], rDZ[7], rDZ[8], DY[2], DY[5], DY[8] );

			//if ( index == INDEX( _nx/2, _ny / 2, _nz - 1 ) )
			//	printf("_rDZ_DX = %3.10e\n", _rDZ_DX.M33[indexOnSurf]);
		}
		else
		{
			pos = INDEX( i, j, 2 * ( _nz - 1 ) - k );
			Jac[index] = Jac[pos];
			con.xi_x[index] = con.xi_x[pos]; con.xi_y[index] = con.xi_y[pos]; con.xi_z[index] = con.xi_z[pos];
			con.et_x[index] = con.et_x[pos]; con.et_y[index] = con.et_y[pos]; con.et_z[index] = con.et_z[pos];
			con.zt_x[index] = con.zt_x[pos]; con.zt_y[index] = con.zt_y[pos]; con.zt_z[index] = con.zt_z[pos];

		}

	END_CALCULATE3D( )

}




void mpiSendRecvJac( GRID grid, MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, float * jac, SEND_RECV_DATA sr );

#ifdef FREE_SURFACE
void solveContravariantJac( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, SEND_RECV_DATA sr, CONTRAVARIANT con, COORD coord, float * Jac, 
MEDIUM medium, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY )
#else
void solveContravariantJac( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, SEND_RECV_DATA sr, CONTRAVARIANT con, COORD coord, float * Jac )
#endif
{
	float DH = grid.DH;
	float rDH = 1.0 / DH;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;




#ifdef GPU_CUDA
	int nx = grid.nx;
	int ny = grid.ny;
	int nz = grid.nz;

	dim3 threads( 32, 4, 4);
	dim3 blocks;
	blocks.x = ( nx + threads.x - 1 ) / threads.x;
	blocks.y = ( ny + threads.y - 1 ) / threads.y;
	blocks.z = ( nz + threads.z - 1 ) / threads.z;

	solve_con_jac
	<<< blocks, threads >>>
	( con, coord, Jac, _nx_, _ny_, _nz_, rDH );

#else
	solve_con_jac( con, coord, Jac, _nx_, _ny_, _nz_, rDH );
#endif
	

#ifdef FREE_SURFACE
	
#ifdef GPU_CUDA
	dim3 threadSurf( 32, 16, 1);
	dim3 blockSurf;
	blockSurf.x = ( nx + threadSurf.x - 1 ) / threadSurf.x;
	blockSurf.y = ( ny + threadSurf.y - 1 ) / threadSurf.y;
	blockSurf.z = HALO + 1;

	solve_coordinate_on_free_surface
	<<< blockSurf, threadSurf >>>
	(	con, coord, Jac, medium, 
		_rDZ_DX, _rDZ_DY, 
		_nx_, _ny_, _nz_); 
#else
	solve_coordinate_on_free_surface
	( con, coord, Jac, medium, 
	  _rDZ_DX, _rDZ_DY, 
	  _nx_, _ny_, _nz_); 
	
#endif //GPU_CUDA

#endif //FREE_SURFACE

	WAVE TW;

	TW.Vx  = con.xi_x; 
	TW.Vy  = con.xi_y; 
	TW.Vz  = con.xi_z; 
	TW.Txx = con.et_x;
	TW.Tyy = con.et_y;
	TW.Tzz = con.et_z;
	TW.Txy = con.zt_x;
	TW.Txz = con.zt_y;
	TW.Tyz = con.zt_z;

	MPI_Barrier( comm_cart );
	mpiSendRecv( grid, comm_cart, mpiNeighbor, TW, sr );
	mpiSendRecvJac( grid, comm_cart, mpiNeighbor, Jac, sr );

}


__GLOBAL__
void matrix_float2FLOAT( CONTRAVARIANT_FLOAT con_fp, CONTRAVARIANT con, FLOAT * Jac_fp, float * Jac, int num )
{
#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif

	CALCULATE1D( i, 0, num )	
		con_fp.xi_x[i] = con.xi_x[i]; 
		con_fp.xi_y[i] = con.xi_y[i]; 
		con_fp.xi_z[i] = con.xi_z[i]; 
		con_fp.et_x[i] = con.et_x[i];
		con_fp.et_y[i] = con.et_y[i];
		con_fp.et_z[i] = con.et_z[i];
		con_fp.zt_x[i] = con.zt_x[i];
		con_fp.zt_y[i] = con.zt_y[i];
		con_fp.zt_z[i] = con.zt_z[i];
		
		Jac_fp[i] = Jac[i];
		
	END_CALCULATE1D( )
}


void matrixfloat2FLOAT( GRID grid, CONTRAVARIANT_FLOAT con_fp, CONTRAVARIANT con, FLOAT * Jac_fp, float * Jac )
{


	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

	long long num = _nx_ * _ny_ * _nz_;

#ifdef GPU_CUDA

	dim3 threadX;
	dim3 blockX;
	threadX.x = 512;
	threadX.y = 1;
	threadX.z = 1;
	blockX.x = ( num + threadX.x - 1 ) / threadX.x;
	blockX.y = 1;
	blockX.z = 1;

	matrix_float2FLOAT
	<<< blockX, threadX >>>
	( con_fp, con, Jac_fp, Jac, num );
#endif

}






/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:station.cpp
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-11-22
*   Discription:
*
================================================================*/
#include "header.h"
void allocStation( STATION * station, int stationNum, int NT )
{
	int sizeIdx = sizeof( int ) * stationNum * 3;

	int * pIndex = NULL;

	CHECK( Malloc( ( void ** )&pIndex, sizeIdx ) );
	CHECK( Memset( pIndex, 0, sizeIdx ) );

	station->X = pIndex;
	station->Y = pIndex + stationNum;
	station->Z = pIndex + stationNum * 2;

	
	long long sizeWave = sizeof( float ) * NT * stationNum * WAVESIZE;
	float * pWave = NULL;
	
	CHECK( Malloc( ( void ** )&pWave, sizeWave ) );
	CHECK( Memset( pWave, 0, sizeWave ) );

	station->wave.Vx  = pWave + NT * stationNum * 0; 
	station->wave.Vy  = pWave + NT * stationNum * 1; 
	station->wave.Vz  = pWave + NT * stationNum * 2; 
	station->wave.Txx = pWave + NT * stationNum * 3;
	station->wave.Tyy = pWave + NT * stationNum * 4;
	station->wave.Tzz = pWave + NT * stationNum * 5;
	station->wave.Txy = pWave + NT * stationNum * 6;
	station->wave.Txz = pWave + NT * stationNum * 7;
	station->wave.Tyz = pWave + NT * stationNum * 8;

}


void allocStation_cpu( STATION * station, int stationNum, int NT )
{
	int sizeIdx = sizeof( int ) * stationNum * 3;

		
	int * pIndex = NULL;

	pIndex = ( int * ) malloc( sizeIdx );
	memset( pIndex, 0, sizeIdx );

	station->X = pIndex;
	station->Y = pIndex + stationNum;
	station->Z = pIndex + stationNum * 2;


	//printf( "alloc:X = %p\n", station->X  );
	//printf( "alloc:Y = %p\n", station->Y  );
	//printf( "alloc:Z = %p\n", station->Z  );
	
	long long sizeWave = sizeof( float ) * NT * stationNum * WAVESIZE;
	float * pWave = NULL;
	
	pWave = ( float * ) malloc( sizeWave );
	memset( pWave, 0, sizeWave );

	station->wave.Vx  = pWave + NT * stationNum * 0; 
	station->wave.Vy  = pWave + NT * stationNum * 1; 
	station->wave.Vz  = pWave + NT * stationNum * 2; 
	station->wave.Txx = pWave + NT * stationNum * 3;
	station->wave.Tyy = pWave + NT * stationNum * 4;
	station->wave.Tzz = pWave + NT * stationNum * 5;
	station->wave.Txy = pWave + NT * stationNum * 6;
	station->wave.Txz = pWave + NT * stationNum * 7;
	station->wave.Tyz = pWave + NT * stationNum * 8;

}

void freeStation( STATION station )
{
	Free( station.X );
	Free( station.wave.Vx );
}

void freeStation_cpu( STATION station )
{
	free( station.X );
	free( station.wave.Vx );
}

int readStationIndex( GRID grid ) 
{
	char jsonFile[1024] = { 0 };
	strcpy( jsonFile, "station.json" );
	FILE * fp;
	fp = fopen( jsonFile, "r" );

	if ( NULL == fp )
	{
		printf( "There is not %s file!\n", jsonFile );
		return 0;
	}
	
	fseek( fp, 0, SEEK_END );
	int len = ftell( fp );
	
	fseek( fp, 0, SEEK_SET );

	char * jsonStr = ( char * ) malloc( len * sizeof( char ) );

	if ( NULL == jsonStr )
	{
		printf( "Can't allocate json string memory\n" );
		return 0;
	}

	fread( jsonStr, sizeof( char ), len, fp );
	
	//printf( "%s\n", jsonStr );
	cJSON * object;
	cJSON * objArray;

	object = cJSON_Parse( jsonStr );
	if ( NULL == object )
	{
		printf( "Can't parse json file!\n");
		//exit( 1 );	
		return 0;
	}

	fclose( fp );



	int stationCnt= 0;

	if (objArray = cJSON_GetObjectItem(object, "station(point)"))
	{
		stationCnt = cJSON_GetArraySize( objArray );
	}
	
	
	cJSON *stationObj, *stationItem;


	int i, j, stationNum;
	int X, Y, Z, thisX, thisY, thisZ;
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;


	stationNum = 0;
	for ( i = 0; i < stationCnt; i ++  )
	{
		stationObj = cJSON_GetArrayItem( objArray, i );


		int a = cJSON_GetArraySize( stationObj );
		if ( a != 3 )
		{
			printf( "In file %s, the coodinate index don't equal to 3. However, it equals to %d\n", jsonFile, a );
			return 0;
		}
	
		stationItem = cJSON_GetArrayItem( stationObj, 0 );
		X = stationItem->valueint;
		thisX = X - frontNX + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 1 );
		Y = stationItem->valueint;
		thisY = Y - frontNY + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 2 );
		Z = stationItem->valueint;
		thisZ = Z - frontNZ + HALO;
			
		if ( thisX >= HALO && thisX < _nx &&
			 thisY >= HALO && thisY < _ny &&
			 thisZ >= HALO && thisZ < _nz )
		{
			stationNum ++;
		}

	}
	//printf( "stationNum = %d\n", stationNum );
	
	return stationNum;


}


void initStationIndex( GRID grid, STATION station ) 
{
	char jsonFile[1024] = { 0 };
	strcpy( jsonFile, "station.json" );
	FILE * fp;
	fp = fopen( jsonFile, "r" );

	if ( NULL == fp )
	{
		printf( "There is not %s file!\n", jsonFile );
		return;
	}
	
	fseek( fp, 0, SEEK_END );
	int len = ftell( fp );
	
	fseek( fp, 0, SEEK_SET );

	char * jsonStr = ( char * ) malloc( len * sizeof( char ) );

	if ( NULL == jsonStr )
	{
		printf( "Can't allocate json string memory\n" );
		return;
	}

	fread( jsonStr, sizeof( char ), len, fp );
	
	//printf( "%s\n", jsonStr );
	cJSON * object;
	cJSON * objArray;

	object = cJSON_Parse( jsonStr );
	if ( NULL == object )
	{
		printf( "Can't parse json file!\n");
		//exit( 1 );	
		return;
	}

	fclose( fp );

	int stationCnt= 0;

	if (objArray = cJSON_GetObjectItem(object, "station(point)"))
	{
		stationCnt = cJSON_GetArraySize( objArray );
		//printf( "stationCnt = %d\n", stationCnt );
	}

	cJSON *stationObj, *stationItem;
	int i, j;
	int thisX, thisY, thisZ;
	
	int frontNX = grid.frontNX;
	int frontNY = grid.frontNY;
	int frontNZ = grid.frontNZ;

	int _nx = grid._nx;
	int _ny = grid._ny;
	int _nz = grid._nz;
	
	int X, Y, Z;



	//printf( "X = %p\n", station.X  );
	//printf( "Y = %p\n", station.Y  );
	//printf( "Z = %p\n", station.Z  );

	int stationIdx = 0;
	for ( i = 0; i < stationCnt; i ++  )
	{
		stationObj = cJSON_GetArrayItem( objArray, i );

		int a = cJSON_GetArraySize( stationObj );
		if ( a != 3 )
		{
			printf( "In file %s, the coodinate index don't equal to 3. However, it equals to %d\n", jsonFile, a );
			return;
		}
	
		stationItem = cJSON_GetArrayItem( stationObj, 0 );
		X = stationItem->valueint;
		thisX = X - frontNX + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 1 );
		Y = stationItem->valueint;
		thisY = Y - frontNY + HALO;

		stationItem = cJSON_GetArrayItem( stationObj, 2 );
		Z = stationItem->valueint;
		thisZ = Z - frontNZ + HALO;
			
		if ( thisX >= HALO && thisX < _nx &&
			 thisY >= HALO && thisY < _ny &&
			 thisZ >= HALO && thisZ < _nz )
		{
			//printf( "X = %p\n", station.X  );
			//printf( "Y = %p\n", station.Y  );
			//printf( "Z = %p\n", station.Z  );
		

			station.X[stationIdx] = thisX;
			station.Y[stationIdx] = thisY;
			station.Z[stationIdx] = thisZ;
		
			stationIdx ++;
		}

	}



}


void stationCPU2GPU( STATION station, STATION station_cpu, int stationNum )
{
	int size = sizeof( int ) * stationNum * 3;
	//int i =0;
	//for ( i = 0; i < stationNum; i ++ )
	//{
	//	printf( "X = %d, Y = %d, Z = %d\n", station_cpu.X[i], station_cpu.Y[i], station_cpu.Z[i] );
	//}

	//printf( "=============size = %d =========\n", size );
	CHECK( Memcpy( station.X, station_cpu.X, size, cudaMemcpyHostToDevice ) );
}


__GLOBAL__
void storage_station( int stationNum, STATION station, WAVE W, int _nx_, int _ny_, int _nz_, int NT, int it )
{

#ifdef GPU_CUDA
	int i = threadIdx.x + blockIdx.x * blockDim.x;
#else
	int i = 0;
#endif
	long long index = 0, pos = 0;
	int X, Y, Z;

	float C1 = 1.0 / Cv;
	float C2 = 1.0 / Cs;

	CALCULATE1D( i, 0, stationNum ) 
		X = station.X[i];
		Y = station.Y[i];
		Z = station.Z[i];
		
		//printf( "X = %d, Y = %d, Z = %d\n", station.X[i], station.Y[i], station.Z[i] );

		//if ( it  % 100 )
		//{
		//	printf( "X = %d, Y = %d, Z = %d\n", X, Y, Z );
		//}

		index = INDEX( X, Y, Z );

		pos = it + i * NT;

		station.wave.Vx [pos] = W.Vx [index] * C1;
		station.wave.Vy [pos] = W.Vy [index] * C1;
		station.wave.Vz [pos] = W.Vz [index] * C1;
		station.wave.Txx[pos] = W.Txx[index] * C2;
		station.wave.Tyy[pos] = W.Tyy[index] * C2;
		station.wave.Tzz[pos] = W.Tzz[index] * C2;
		station.wave.Txy[pos] = W.Txy[index] * C2;
		station.wave.Txz[pos] = W.Txz[index] * C2;
		station.wave.Tyz[pos] = W.Tyz[index] * C2;
	
	END_CALCULATE1D( )

}

void storageStation( GRID grid, int NT, int stationNum, STATION station, WAVE W, int it )
{
	long long num = stationNum;

	int _nx_ = grid._nx_;
	int _ny_ = grid._ny_;
	int _nz_ = grid._nz_;

#ifdef GPU_CUDA
	dim3 threads( 8, 8, 8 );
	dim3 blocks;
	blocks.x = ( num + threads.x - 1 ) / threads.x;
	blocks.y = 1;
	blocks.z = 1;
	storage_station
	<<< blocks, threads >>>
	( stationNum, station, W, _nx_, _ny_, _nz_, NT, it );
	CHECK( cudaDeviceSynchronize( ) );
#else
	storage_station
	( stationNum, station, W, _nx_, _ny_, _nz_, NT, it );
#endif
}

void stationGPU2CPU( STATION station, STATION station_cpu, int stationNum, int NT )
{
	long long sizeWave = sizeof( float ) * NT * stationNum * WAVESIZE;
	CHECK( Memcpy( station_cpu.wave.Vx, station.wave.Vx, sizeWave, cudaMemcpyDeviceToHost ) );
}

void write( PARAMS params, GRID grid, MPI_COORD thisMPICoord, STATION station, int NT, int stationNum )
{	
	FILE * fp;
	char fileName[256];
	sprintf( fileName, "%s/station_mpi_%d_%d_%d.bin", params.OUT, thisMPICoord.X, thisMPICoord.Y, thisMPICoord.Z );


	int i = 0;

	for ( i = 0; i < stationNum; i ++ )
	{
		station.X[i] = grid.frontNX + station.X[i] - HALO;
		station.Y[i] = grid.frontNY + station.Y[i] - HALO;
		station.Z[i] = grid.frontNZ + station.Z[i] - HALO;
	}

	fp = fopen( fileName, "wb" ); 

	fwrite( station.X, sizeof( int ), stationNum * 3, fp );
	fwrite( station.wave.Vx, sizeof( float ), NT * stationNum * WAVESIZE, fp );
	
	fclose( fp );

}




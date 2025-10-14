#include "knn.cuh"

#include <utility>
#include "Vortex2D.h"
#include <omp.h>
#include <algorithm>
#include <iostream>

#include "wrapper.h"




struct prDoubleInt
{
	double first;
	size_t second;
};

class VectorsForKnn
{
public:
	VectorsForKnn();
	~VectorsForKnn();

	//указатель на вихри
	Vortex2D* vtxPtr; 
		
	//мортоновские коды вихрей
	int* mcdataPtr;
	int* mcdata_unsortedPtr;

	//пор€док сортировки
	int* indexPtr;
	int* index_unsortedPtr;

	prDoubleInt* initdistPtr;       // указатели на структуру данных соседей = { рассто€ние; до какого элемента (в отсортированном массиве) }
	prDoubleInt* initdistPtrSdvig;  // то же, но после сдвига мортоновских кодов, потом мерджитс€ initdistPtr
	prDoubleInt* initdistPtrUpdate; // то же, имеет смысл временной переменной дл€ мерджинга
	
	int* loc;
	int* offset;
	int* counter;
	int* counterskan;

	prDoubleInt* dstKeys1Ptr; //временна€ переменна€ дл€ сортировки

	void* sortBuffer;
	int sortBufferSizeInBytes;

	void ResizeAll(int k, int nvtx);

private:
	int reservedVtx;

};


VectorsForKnn::VectorsForKnn() : reservedVtx(0), sortBufferSizeInBytes(0) {};

VectorsForKnn::~VectorsForKnn() {
	cudaFree(initdistPtr);
	cudaFree(initdistPtrSdvig);
	cudaFree(initdistPtrUpdate);
	cudaFree(loc);
	cudaFree(offset);
	cudaFree(counter);
	cudaFree(counterskan);

	cudaFree(vtxPtr);

	cudaFree(mcdataPtr);
	cudaFree(mcdata_unsortedPtr);

	cudaFree(indexPtr);
	cudaFree(index_unsortedPtr);

	cudaFree(dstKeys1Ptr);

	cudaFree(sortBuffer);
};

void VectorsForKnn::ResizeAll(int k, int nvtx)
{
	if ((nvtx > 0) && (nvtx > reservedVtx))
	{
		if (reservedVtx > 0)
		{
			cudaFree(vtxPtr);
			cudaFree(mcdataPtr);
			cudaFree(mcdata_unsortedPtr);
			cudaFree(indexPtr);
			cudaFree(index_unsortedPtr);
			cudaFree(dstKeys1Ptr);
			cudaFree(initdistPtr);
			cudaFree(initdistPtrSdvig);
			cudaFree(initdistPtrUpdate);
			cudaFree(loc);
			cudaFree(offset);
			cudaFree(counter);
			cudaFree(counterskan);
			cudaFree(sortBuffer);
		}

		reservedVtx = nvtx * 2;
		std::cout << "knn: reservedVtx = " << reservedVtx << ", nvtx = " << nvtx << std::endl;

		cudaMalloc((void**)&vtxPtr, reservedVtx * sizeof(Vortex2D));		
		
		cudaMalloc((void**)&mcdataPtr, reservedVtx * sizeof(int));
		cudaMalloc((void**)&mcdata_unsortedPtr, reservedVtx * sizeof(int));		
		
		cudaMalloc((void**)&indexPtr, reservedVtx * sizeof(int));
		cudaMalloc((void**)&index_unsortedPtr, reservedVtx * sizeof(int));

		cudaMalloc((void**)&dstKeys1Ptr, 2 * k * reservedVtx * sizeof(prDoubleInt));
		cudaMalloc((void**)&initdistPtr, reservedVtx * 2 * k * sizeof(prDoubleInt));
		cudaMalloc((void**)&initdistPtrSdvig, reservedVtx * 2 * k * sizeof(prDoubleInt));
		cudaMalloc((void**)&initdistPtrUpdate, reservedVtx * 2 * k * sizeof(prDoubleInt));

		cudaMalloc((void**)&loc, reservedVtx * 2 * k * sizeof(int));
		cudaMalloc((void**)&offset, reservedVtx * 2 * k * sizeof(int));
		cudaMalloc((void**)&counter, reservedVtx * 2 * k * sizeof(int));
		cudaMalloc((void**)&counterskan, reservedVtx * 2 * k * sizeof(int));

		cudaMalloc((void**)&sortBuffer, reservedVtx * 10 * sizeof(char));
		sortBufferSizeInBytes = reservedVtx * 10;
	}
}


const int codeLength = 15;
const int twoPowCodeLength = (1 << codeLength);

//"–азрежение" двоичного представлени€ беззнакового целого, вставл€€ по одному нулику между всеми битами
__host__ __device__ unsigned int ExpandBitscuda(unsigned int v)
{
	// вставит 1 нуль
	v = (v | (v << 8)) & 0x00FF00FF;      //  00000000`00000000`abcdefgh`ijklmnop 
	//                                      | 00000000`abcdefgh`ijklmnop`00000000
	//                                      = 00000000`abcdefgh`XXXXXXXX`ijklmnop
	//                                      & 00000000`11111111`00000000`11111111
	//                                      = 00000000`abcdefgh`00000000`ijklmnop

	v = (v | (v << 4)) & 0x0F0F0F0F;      //  00000000`abcdefgh`00000000`ijklmnop 
	//                                      | 0000abcd`efgh0000`0000ijkl`mnop0000
	//                                      = 0000abcd`XXXXefgh`0000ijkl`XXXXmnop
	//                                      & 00001111`00001111`00001111`00001111
	//                                      = 0000abcd`0000efgh`0000ijkl`0000mnop

	v = (v | (v << 2)) & 0x33333333;      //  0000abcd`0000efgh`0000ijkl`0000mnop 
	//                                      | 00abcd00`00efgh00`00ijkl00`00mnop00
	//                                      = 00abXXcd`00efXXgh`00ijXXkl`00mnXXop
	//                                      & 00110011`00110011`00110011`00110011
	//                                      = 00ab00cd`00ef00gh`00ij00kl`00mn00op

	v = (v | (v << 1)) & 0x55555555;      //  00ab00cd`00ef00gh`00ij00kl`00mn00op 
	//                                      | 0ab00cd0`0ef00gh0`0ij00kl0`0mn00op0
	//                                      = 0aXb0cXd`0eXf0gXh`0iXj0kXl`0mXn0oXp
	//                                      & 01010101`01010101`01010101`01010101
	//                                      = 0a0b0c0d`0e0f0g0h`0i0j0k0l`0m0n0o0p
	return v;
}


//ћортоновский код дл€ пары из чисел типа double
//»сходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
__host__ __device__ int Morton2Dcuda(const Point2D& r)
{
	const Point2D& rscale = twoPowCodeLength * r;
	const unsigned int& xx = ExpandBitscuda((unsigned int)(rscale[0]));
	const unsigned int& yy = ExpandBitscuda((unsigned int)(rscale[1]));
	return (int)(yy | (xx << 1));
}


template <typename IT>
__host__ __device__ size_t BinSearchcuda(
	IT beg, IT end,
	double x, int low, int high)
{
	int mid = -1;

	if (x > (beg + high)->first)
		return high + 1;

	while (low <= high) {
		mid = (low + high) / 2;
		if ((beg + mid)->first == x)
			return mid + 1; //выход из цикла
		if ((beg + mid)->first < x)
			low = mid + 1;
		else
			high = mid - 1;
	}
	return mid;
}


template <typename IT, typename IT2>
__host__ __device__ void newSortcuda(IT beg, IT end, IT2 dstKeys) {
	size_t k = end - beg;
	size_t cnt = 0;

	for (size_t i = 0; i < k; i++) {
		double elem = (beg + i)->first;

		cnt = 0;
		for (size_t j = 0; j < k; j++) {
			cnt += ((beg + j)->first < elem);
		}
		//if (!AssumeDistinct) 
		for (size_t j = 0; j < i; j++)
			cnt += ((beg + j)->first == elem);
		*(dstKeys + cnt) = *(beg + i);
	}

	for (int i = 0; i < k; ++i)
	{
		*(beg + i) = *(dstKeys + i);
	}
}

template <typename IT, typename IT2, typename IT3>
void __host__ __device__ newMergecuda(
	int iii,
	size_t k,
	IT beg, IT end,	
	IT3 candidateNN,		 
	IT2 loc,  //2k
	IT2 counter,//2k
	IT2 offset,//2k
	IT2 counterScan,//2k
	IT3 updateNN//k
)
{
	for (size_t j = 0; j < 2 * k; ++j)
		*(loc + j) = BinSearchcuda(beg, end, (*(candidateNN+j)).first, 0, (int)k - 1);

	for (size_t j = 0; j < 2 * k; ++j)
		*(counter + j) = j & 1;

	for (size_t j = 0; j < 2 * k; ++j)
	{		
		if (( *(loc+j) > 0) && ( (*(loc+j) == k) || ((*(candidateNN + j)).second == (beg + *(loc + j) - 1)->second)))
			*(offset+j) = k+1;
		else
			*(offset+j) = (*(counter + 2 * *(loc+j)))++;
	}

	*(counterScan+0) = 0;
	for (size_t j = 1; j < 2 * k; ++j)
		*(counterScan+j) = *(counterScan + j-1) + *(counter + j - 1);


	size_t index;
	for (size_t j = 0; j < k; ++j)
	{
		index = *(counterScan + 2 * j + 1);
		if (index < k)
			*(updateNN+index) = *(beg + j);
	}

	for (size_t j = 0; j < 2 * k; ++j)
	{
		if (2 * *(loc+j) < 2 * k)
		{
			index = counterScan[2 * *(loc + j)] + *(offset+j);
			if (index < k)
				*(updateNN + index) = candidateNN[j];
		}
	}

	for (int s = 0; s < k; ++s)
	{
		(*(beg + s)).first = (*(updateNN + s)).first;
		(*(beg + s)).second = (*(updateNN + s)).second;
	}
}

#include "cuSort.cuh"


__device__ void calcCheck(bool ch, const Vortex2D& vtxi, const Vortex2D& vtxk, double maxG, double cSP, double cRBP, double epsCol, int type, double& d2, int& result1, int& result2)
{
	int flagExit = false;
	int check = false;

	//линейное увеличение радиуса коллапса				    
	double mnog = std::max(1.0, /* 2.0 * */ (vtxi.r()[0] - cRBP) / cSP);
	double r2test = (epsCol * mnog) * (epsCol * mnog);

	if (type == 1)
		r2test *= 4.0; //”величение радиуса коллапса в 2 раза дл€ коллапса вихрей разных знаков		

	d2 = (vtxi.r() - vtxk.r()).length2();

	if (d2 > 1.0 * r2test)
	{
		flagExit = true;
		result1 = flagExit;
		result2 = false;
		return;
	}

	const double gi = vtxi.g();
	const double gk = vtxk.g();

	if (d2 < r2test)
	{
		switch (type)
		{
		case 0:		
			check = (fabs(gi * gk) != 0.0) && (fabs(gi + gk) < (mnog * mnog) * maxG);
			break;
		case 1:
			check = (gi * gk < 0.0);
			break;
		case 2:
			check = (gi * gk > 0.0) && (fabs(gi + gk) < (mnog * mnog) * maxG);
			break;
		}
	}//if r2 < r2_test

	result1 = flagExit;
	result2 = check;
	return;
}




__global__ void findNeib(
	prDoubleInt* initdistPtr,
	prDoubleInt* initdistPtrSdvig,
	prDoubleInt* initdistPtrUpdate,
	int* loc, 
	int* offset, 
	int* counter, 
	int* counterskan,	
	int* index_dev_raw, 
	Vortex2D* devVtxPtrV2D, 
	int k, 
	int n, 
	int sdvig, 
	double maxG, 
	double cSP,
	double cRBP, 
	double epsCol, 
	int type,
	prDoubleInt* dstKeys1Ptr
	)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	prDoubleInt* curPosition = nullptr;

	int flagExit, check;
	double d2;

	int cntr = 0;
	int search = i - 1;

	if (i < n)
	{
		const Vortex2D& vtxi = devVtxPtrV2D[index_dev_raw[i]];

		while ((cntr < k) && (search >= 0))
		{
			const Vortex2D& vtxk = devVtxPtrV2D[index_dev_raw[search]];
			calcCheck(false, vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2, flagExit, check);
			if (flagExit)
				break;

			if ((index_dev_raw[search] > index_dev_raw[i]) && check)
			{
				//__syncthreads();
				curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (index_dev_raw[i] * 2 * k) + cntr;
				curPosition->first = d2;
				curPosition->second = index_dev_raw[search];
				++cntr;
			}
			--search;
		}

		__threadfence();
		__syncthreads();


		for (int w = cntr; w < k; ++w)			
		{
			curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (index_dev_raw[i] * 2 * k) + cntr;					
			curPosition->first = 100000000.0;
			curPosition->second = 0;
			++cntr;			
		}
		__threadfence();
		__syncthreads();

		search = i + 1;

		while ((cntr < 2 * k) && (search < n))
		{
			const Vortex2D& vtxk = devVtxPtrV2D[index_dev_raw[search]];
			calcCheck(false, vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2, flagExit, check);
			if (flagExit)
				break;

			if ((index_dev_raw[search] > index_dev_raw[i]) && check)
			{
				curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (index_dev_raw[i] * 2 * k) + cntr;
				curPosition->first = d2;
				curPosition->second = index_dev_raw[search];
				++cntr;
			}
			++search;
		}
		__threadfence();
		__syncthreads();

		for (int w = cntr; w < 2 * k; ++w)
		{
			curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (index_dev_raw[i] * 2 * k) + cntr;
			curPosition->first = 100000000.0;
			curPosition->second = 0;
			++cntr;			
		}
		__threadfence();
		__syncthreads();

		if (sdvig == 0)					
			newSortcuda(initdistPtr + (index_dev_raw[i] * 2 * k), initdistPtr + (index_dev_raw[i] + 1) * 2 * k, dstKeys1Ptr + (index_dev_raw[i] * 2 * k));
		else
		{
			newMergecuda(index_dev_raw[i], k,
				initdistPtr + (index_dev_raw[i] * 2 * k), initdistPtr + (index_dev_raw[i]) * 2 * k + k,
				initdistPtrSdvig + (index_dev_raw[i] * 2 * k),
				loc + (index_dev_raw[i] * 2 * k),
				counter + (index_dev_raw[i] * 2 * k),
				offset + (index_dev_raw[i] * 2 * k),
				counterskan + (index_dev_raw[i] * 2 * k),
				initdistPtrUpdate + (index_dev_raw[i] * 2 * k)
			);
			newSortcuda(initdistPtr + (index_dev_raw[i] * 2 * k), initdistPtr + (index_dev_raw[i] + 1) * 2 * k, dstKeys1Ptr + (index_dev_raw[i] * 2 * k));
		}
		__threadfence();
		__syncthreads();
	}
};


__global__ void Fill123(int* ptr, int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n)	
		ptr[i] = i;
};



__global__ void CalcMortonCodes(Vortex2D* vtxPtr, int n, Point2D LL, double scale, int sdvig, Point2D shift005, int* mcdata_unsortedPtr)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n)
	{
		Point2D sh = ((vtxPtr[i].r() - LL) * (0.75 / scale) + sdvig * shift005);
		mcdata_unsortedPtr[i] = Morton2Dcuda(sh);
	};
}

__global__ void reordNeib(int n, int k,
	prDoubleInt* initdistPtr,
	prDoubleInt* initdistPtrUpdate)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n)
	{
		for (int q = 0; q < k; ++q)
			initdistPtrUpdate[i * k + q] = initdistPtr[i * (2 * k) + q];

	}//if n
}


void kNNcuda(const std::vector<Vortex2D>& vtx,
	const size_t k,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int type
)
{
	//double t1 = -omp_get_wtime();
	const size_t n = vtx.size();

	if (vecForKnn == nullptr)
		vecForKnn = new VectorsForKnn;

	vecForKnn->ResizeAll((int)k, (int)n);
	cudaMemcpy(vecForKnn->vtxPtr, vtx.data(), n * sizeof(Vortex2D), cudaMemcpyHostToDevice);

	//ѕоиск габарита
	BHcu::CudaCalcGab gab;
	Point2D minr, maxr;
	gab.calc((int)n, vecForKnn->vtxPtr);

	cudaMemcpy(&minr, gab.minpt, sizeof(Point2D), cudaMemcpyDeviceToHost);
	cudaMemcpy(&maxr, gab.maxpt, sizeof(Point2D), cudaMemcpyDeviceToHost);

	double scale = std::max(maxr[0] - minr[0], maxr[1] - minr[1]);
	
	const int nSdvig = 5;
	Point2D shift005;
	shift005[0] = shift005[1] = 0.05;

	dim3 blockA, threadA = 32;
	blockA.x = ((int)n + 31) / 32;

	double tA = omp_get_wtime();

	for (size_t sdvig = 0; sdvig < nSdvig; ++sdvig)
	{
		//масштабирование координат вихрей в [0;0.75)^2, поиск их мортоновских кодов 	

		//—двигаем все точки и считаем новые коды
		CalcMortonCodes<<<blockA, threadA>>> (vecForKnn->vtxPtr, (int)n, minr, scale, (int)sdvig, shift005, vecForKnn->mcdata_unsortedPtr);

		Fill123<<<blockA, threadA>>>(vecForKnn->index_unsortedPtr, (int)n);

		BHcu::RadixSortFromCUBReservedMem(
			vecForKnn->mcdata_unsortedPtr,
			vecForKnn->mcdataPtr,
			vecForKnn->index_unsortedPtr,
			vecForKnn->indexPtr,
			(int)(n), 0, 2 * codeLength,
			vecForKnn->sortBuffer,
			vecForKnn->sortBufferSizeInBytes
			);

		findNeib<<<blockA, threadA>>>(
			vecForKnn->initdistPtr,
			vecForKnn->initdistPtrSdvig,
			vecForKnn->initdistPtrUpdate,
			vecForKnn->loc,
			vecForKnn->offset,
			vecForKnn->counter,
			vecForKnn->counterskan,			
			vecForKnn->indexPtr,
			vecForKnn->vtxPtr,			
			(int)k,
			(int)n,
			(int)sdvig,
			maxG,
			cSP,
			cRBP,
			epsCol,
			type,
			vecForKnn->dstKeys1Ptr);
	}//for sdvig

	double tB = omp_get_wtime();
	//std::cout << "tB-tA = " << tB-tA << std::endl;

	reordNeib << <blockA, threadA >> > ((int)n, (int)k, vecForKnn->initdistPtr, vecForKnn->initdistPtrUpdate);

	double tMem = -omp_get_wtime();
	cudaMemcpy(initdist.data(), (char*)vecForKnn->initdistPtrUpdate, n * k * sizeof(prDoubleInt), cudaMemcpyDeviceToHost);

	tMem +=omp_get_wtime(); 
	//std::cout << "Tmem = " << tMem << std::endl;
}

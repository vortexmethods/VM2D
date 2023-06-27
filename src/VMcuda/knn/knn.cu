#include "knn.cuh"

#include <utility>
#include "Vortex2D.h"
#include <omp.h>
#include <algorithm>


#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include <thrust/partition.h>

class VectorsForKnn
{
public:
	VectorsForKnn();
	~VectorsForKnn();

	thrust::device_vector<Vortex2D> vtx_dev;
	//для каждого запроса запоминаем индекс iq -- индекс ближайшего левого элемента данных
	//отсортиванный массив mortonCodeBoth разделяем на mcquery (только запросы) и mcdata (только данные)
	thrust::device_vector<TParticleCode> mcdatadev, mcquerydev;
	thrust::device_vector<unsigned int> iqdev;

	//поиск для каждого запроса k точек данных слева и k точек данных справа (при нехватке данных с одной стороны, добавляем с другой)
	//для каждого запроса находим индекс элемента данных, с которого будут отсчитываться 2*k соседей
	thrust::device_vector<size_t> initneigdev;

	//выбор из 2*k кандидатов k соседей (выбор по расстоянию от данных до запроса)
	// initdist = {расстояние; до элемента данных}
	// сортировка по расстоянию

	//Генерируем мортоновы коды для точек и запросов (идентичны)
	thrust::device_vector<int> unsorted_mortonCodeBoth_dev;
	thrust::device_vector<int> mortonCodeBoth_dev;

	thrust::device_vector<int> unsorted_indicesBoth_dev;
	thrust::device_vector<int> indicesBoth_dev;

	thrust::device_vector<TParticleCode> mortonCodeIndexBoth_dev;

	thrust::device_vector<int> prfx;
	thrust::device_vector<int> posTrue;
	thrust::device_vector<thrust::pair<double, size_t>> dstKeys1_dev;
	thrust::device_vector<thrust::pair<double, size_t>> dstKeys_dev;

	thrust::device_vector<thrust::pair<double, size_t>> initdist_dev;
	thrust::device_vector<thrust::pair<double, size_t>> dist_dev;

	thrust::device_vector<size_t> loc_dev;
	thrust::device_vector<size_t> counter_dev;
	thrust::device_vector<size_t> offset_dev;
	thrust::device_vector<size_t> counterScan_dev;
	thrust::device_vector<thrust::pair<double, size_t>> updateNN_dev;

	void ResizeAll(int k, int nvtx);

private:
	int reservedVtx;

};


VectorsForKnn::VectorsForKnn() : reservedVtx(0) {};
VectorsForKnn::~VectorsForKnn() {};

void VectorsForKnn::ResizeAll(int k, int nvtx)
{

	if ((reservedVtx == 0) || (nvtx > reservedVtx))
	{
		reservedVtx = nvtx * 2;
		std::cout << "reservedVtx = " << reservedVtx << ", nvtx = " << nvtx << std::endl;

		vtx_dev.resize(reservedVtx);

		//для каждого запроса запоминаем индекс iq -- индекс ближайшего левого элемента данных
		//отсортиванный массив mortonCodeBoth разделяем на mcquery (только запросы) и mcdata (только данные)
		mcdatadev.resize(reservedVtx); 
		mcquerydev.resize(reservedVtx);
		iqdev.resize(reservedVtx);

		//поиск для каждого запроса k точек данных слева и k точек данных справа (при нехватке данных с одной стороны, добавляем с другой)
		//для каждого запроса находим индекс элемента данных, с которого будут отсчитываться 2*k соседей
		initneigdev.resize(reservedVtx, (size_t)(-1));

		//выбор из 2*k кандидатов k соседей (выбор по расстоянию от данных до запроса)
		// initdist = {расстояние; до элемента данных}
		// сортировка по расстоянию

		//Генерируем мортоновы коды для точек и запросов (идентичны)
		unsorted_mortonCodeBoth_dev.resize(2 * reservedVtx);
		mortonCodeBoth_dev.resize(2 * reservedVtx);

		unsorted_indicesBoth_dev.resize(2 * reservedVtx);
		indicesBoth_dev.resize(2 * reservedVtx);

		mortonCodeIndexBoth_dev.resize(2 * reservedVtx);

		prfx.resize(2 * reservedVtx);
		posTrue.resize(reservedVtx);

		dstKeys1_dev.resize(2 * k * reservedVtx, { 0.0, 0 });
		dstKeys_dev.resize(k * reservedVtx, { 0.0, 0 });

		initdist_dev.resize(2 * k * reservedVtx);
		dist_dev.resize(2 * k * reservedVtx);

		loc_dev.resize(2 * k * reservedVtx);
		counter_dev.resize(2 * k * reservedVtx);
		offset_dev.resize(2 * k * reservedVtx);
		counterScan_dev.resize(2 * k * reservedVtx);
		updateNN_dev.resize(k * reservedVtx);
	}
}



const int codeLength = 15;
const int twoPowCodeLength = (1 << codeLength);

//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
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


//Мортоновский код для пары из чисел типа double
//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
__host__ __device__ int Morton2Dcuda(const Point2D& r)
{
	const Point2D& rscale = twoPowCodeLength * r;
	const unsigned int& xx = ExpandBitscuda((unsigned int)(rscale[0]));
	const unsigned int& yy = ExpandBitscuda((unsigned int)(rscale[1]));
	return (int)(yy | (xx << 1));
}


/// Поиск границ габаритного прямоугольника (левый нижний и правый верхний) 
std::pair<Point2D, Point2D> GetGabaritcuda(Vortex2D* vtx_dev_ptr, int nvtx)
{
	thrust::pair < thrust::device_ptr<Vortex2D>, thrust::device_ptr<Vortex2D>> 
		xx = thrust::minmax_element(thrust::device_ptr<Vortex2D>(vtx_dev_ptr), thrust::device_ptr<Vortex2D>(vtx_dev_ptr) + nvtx, []__device__(const Vortex2D & P1, const Vortex2D & P2) { return P1.r()[0] < P2.r()[0]; });
	
	thrust::pair < thrust::device_ptr<Vortex2D>, thrust::device_ptr<Vortex2D>>
		yy = thrust::minmax_element(thrust::device_ptr<Vortex2D>(vtx_dev_ptr), thrust::device_ptr<Vortex2D>(vtx_dev_ptr) + nvtx, []__device__(const Vortex2D & P1, const Vortex2D & P2) { return P1.r()[1] < P2.r()[1]; });

	Vortex2D LLx = *xx.first;
	Vortex2D LLy = *yy.first;

	Vortex2D URx = *xx.second;	
	Vortex2D URy = *yy.second;

	Point2D LL, UR;		
	LL[0] = LLx.r()[0];
	LL[1] = LLy.r()[1];	
	UR[0] = URx.r()[0];
	UR[1] = URy.r()[1];

 	return { LL, UR };
}


template <typename IT>
__host__ __device__ size_t BinSearchcuda(
	IT beg, IT end,
	double x, int low, int high)
{
	int mid = -1;

	//if (x > currentNN[high].first)
	if (x > (beg + high)->first)
		return high;

	while (low <= high) {
		mid = (low + high) / 2;

		//if (currentNN[mid].first == x)
		if ((beg + mid)->first == x)
			return mid + 1; //выход из цикла

		//if (currentNN[mid].first < x)
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
		(*(beg + i)).first = (*(dstKeys + i)).first;
		(*(beg + i)).second = (*(dstKeys + i)).second;
	}
}

template <typename IT, typename IT2, typename IT3>
void __host__ __device__ newMergecuda(
	size_t k,
	IT beg, IT end,
	
	IT3 candidateNN,	 
	 
	IT2 loc,
	IT2 counter,
	IT2 offset,
	IT2 counterScan,

	IT3 updateNN
)
{
	for (size_t j = 0; j < 2 * k; ++j)
		*(loc + j) = BinSearchcuda(beg, end, (*(candidateNN+j)).first, 0, (int)k - 1);

	for (size_t j = 0; j < 2 * k; ++j)
		*(counter + j) = j & 1;

	for (size_t j = 0; j < 2 * k; ++j)
	{		
		if (( *(loc+j) > 0) && ( (*(loc+j) == k) || ((*(candidateNN + j)).first == (beg + *(loc+j) - 1)->first))) //можно через second, проверить!!
			*(offset+j) = k;
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

__host__ void mainCyclecuda(
	const int k,
	const std::vector<Vortex2D>& vtx, 
	const Vortex2D* vtx_dev_raw,
	const std::pair<Point2D, Point2D>& gab,
	const double scale, const int sdvig,

	thrust::device_vector<TParticleCode>& mcdata_dev,
	thrust::device_vector<TParticleCode>& mcquery_dev,
	
	thrust::device_vector<unsigned int>& iq_dev,
	thrust::device_vector<size_t>& initneig_dev,
	double* time,

	thrust::device_vector<int>& unsorted_mortonCodeBoth_dev,
	thrust::device_vector<int>& mortonCodeBoth_dev,

	thrust::device_vector<int>& unsorted_indicesBoth_dev,
	thrust::device_vector<int>& indicesBoth_dev,

	thrust::device_vector<TParticleCode>& mortonCodeIndexBoth_dev,


	thrust::device_vector<int>& prfx,
	thrust::device_vector<int>& posTrue,
	thrust::device_vector<thrust::pair<double, size_t>>& dstKeys1_dev,
	thrust::device_vector<thrust::pair<double, size_t>>& dstKeys_dev,
	thrust::device_vector<thrust::pair<double, size_t>>& initdist_dev,
	thrust::device_vector<thrust::pair<double, size_t>>& dist_dev,
	thrust::device_vector<size_t>& loc_dev,
	thrust::device_vector<size_t>& counter_dev,
	thrust::device_vector<size_t>& offset_dev,
	thrust::device_vector<size_t>& counterScan_dev,
	thrust::device_vector<thrust::pair<double, size_t>>& updateNN_dev
)

{

	thrust::pair<double, size_t>* dstKeys1_dev_raw = thrust::raw_pointer_cast(&dstKeys1_dev[0]);
	thrust::pair<double, size_t>* dstKeys_dev_raw = thrust::raw_pointer_cast(&dstKeys_dev[0]);

	thrust::pair<double, size_t>* initdist_dev_raw = thrust::raw_pointer_cast(&initdist_dev[0]);
	TParticleCode* mcquery_dev_raw = thrust::raw_pointer_cast(&mcquery_dev[0]);
	TParticleCode* mcdata_dev_raw = thrust::raw_pointer_cast(&mcdata_dev[0]);
	size_t* initneig_dev_raw = thrust::raw_pointer_cast(&initneig_dev[0]);



	//time[0] = omp_get_wtime();

	size_t nvtx = vtx.size();

	//масштабирование координат вихрей в [0;0.75)^2, поиск их мортоновских кодов 	
	Point2D shift005;
	shift005[0] = 0.05;
	shift005[1] = 0.05;
	
	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
		unsorted_mortonCodeBoth_dev.begin(),
		[vtx_dev_raw, gab, scale, sdvig, shift005] __device__(const unsigned int idx) {
		int code = Morton2Dcuda((vtx_dev_raw[idx].r() - gab.first) * (0.75 / scale) + sdvig * shift005);
		code <<= 1;	
		return code;
	});

	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
		unsorted_indicesBoth_dev.begin(),
		[] __device__(const unsigned int idx) {		
		return idx;
	});

	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
		unsorted_mortonCodeBoth_dev.begin() + (unsigned int)nvtx,		
		[vtx_dev_raw, gab, scale, sdvig, shift005] __device__(const unsigned int idx) {
		int code = Morton2Dcuda((vtx_dev_raw[idx].r() - gab.first) * (0.75 / scale) + sdvig * shift005);
		code <<= 1;
		code |= 1;		
		return code;
	});

	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
		unsorted_indicesBoth_dev.begin() + (unsigned int)nvtx,
		[] __device__(const unsigned int idx) {
		return idx;
	});

	//time[1] = omp_get_wtime();
	
	//Сортировка мортоновых кодов (используя Cuda)
	BHcu::RadixSortFromCUB(
		thrust::raw_pointer_cast(&unsorted_mortonCodeBoth_dev[0]),
		thrust::raw_pointer_cast(&mortonCodeBoth_dev[0]),
		thrust::raw_pointer_cast(&unsorted_indicesBoth_dev[0]),
		thrust::raw_pointer_cast(&indicesBoth_dev[0]),
		(int)(2 * nvtx), 0, 2 * codeLength);
	
	int* mortonCodeBoth_dev_raw = thrust::raw_pointer_cast(&mortonCodeBoth_dev[0]);
	int* indicesBoth_dev_raw = thrust::raw_pointer_cast(&indicesBoth_dev[0]);

	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)(2 * nvtx)),
		mortonCodeIndexBoth_dev.begin(),
		[mortonCodeBoth_dev_raw, indicesBoth_dev_raw] __device__(const unsigned int idx) {
		TParticleCode code;
		code.key = mortonCodeBoth_dev_raw[idx];
		code.originNumber = indicesBoth_dev_raw[idx];
		return code;
	});

	TParticleCode* mortonCodeIndexBoth_dev_raw = thrust::raw_pointer_cast(&mortonCodeIndexBoth_dev[0]);


	//time[2] = omp_get_wtime();
	//для каждого запроса запоминаем индекс iq -- индекс ближайшего левого элемента данных
    //отсортиванный массив mortonCodeBoth разделяем на mcquery (только запросы) и mcdata (только данные)
	
	thrust::stable_partition_copy(
		mortonCodeIndexBoth_dev.begin(), mortonCodeIndexBoth_dev.begin() + nvtx * 2,
		mcquery_dev.begin(), mcdata_dev.begin(),
		[]__device__(const TParticleCode& p) {
		bool cond = (p.key % 2);
		return (cond);
			// ((unsigned long)p.key & (unsigned long)1); 
	});
	

	thrust::transform(
		mortonCodeIndexBoth_dev.begin(), mortonCodeIndexBoth_dev.begin() + 2 * nvtx,
		prfx.begin(),
		[] __device__(const TParticleCode& p) {
		return !(p.key & 1);		
	});

	thrust::inclusive_scan(prfx.begin(), prfx.begin() + 2 * nvtx, prfx.begin());
	
	thrust::copy_if(
		thrust::make_counting_iterator<int>(0),
		thrust::make_counting_iterator<int>((int)(2 * nvtx)),
		posTrue.begin(),
		[mortonCodeIndexBoth_dev_raw] __device__(const unsigned int idx) {
		return (mortonCodeIndexBoth_dev_raw[idx].key & 1);		
	});

	int* prfx_raw = thrust::raw_pointer_cast(&prfx[0]);
	int* posTrue_raw = thrust::raw_pointer_cast(&posTrue[0]);

	thrust::transform(
		thrust::make_counting_iterator<int>(0),
		thrust::make_counting_iterator<int>((int)(nvtx)),
		iq_dev.begin(),
		[prfx_raw, posTrue_raw]__device__(unsigned int idx)
	{
		return (prfx_raw[posTrue_raw[idx]] == 0) ? 0 : prfx_raw[posTrue_raw[idx]] - 1;
	});

	//time[3] = omp_get_wtime();
	//time[4] = omp_get_wtime();

	//поиск для каждого запроса k точек данных слева и k точек данных справа (при нехватке данных с одной стороны, добавляем с другой)
	//для каждого запроса находим индекс элемента данных, с которого будут отсчитываться 2*k соседей			

	unsigned int* iq_dev_raw = thrust::raw_pointer_cast(&iq_dev[0]);

	thrust::transform(
		thrust::make_counting_iterator<unsigned int>(0),
		thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
		initneig_dev.begin(),
		[iq_dev_raw, k, nvtx]__device__(unsigned int idx)
	{
		size_t pt = iq_dev_raw[idx];
		size_t left = pt - k + 1;
		if (pt < k - 1)
			left = 0;
		if (pt > nvtx - k - 1)
			left = (nvtx - 1) - (2 * k - 1);
		
		return left;
	});


	if (sdvig == 0)
	{								
		thrust::for_each(			
			thrust::make_counting_iterator<unsigned int>(0),
			thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
			[initdist_dev_raw, mcquery_dev_raw, mcdata_dev_raw, vtx_dev_raw, initneig_dev_raw, k, nvtx]
			__device__(unsigned int idx)			
		{
	
			for (size_t j = 0; j < 2 * k; ++j)
			{
				initdist_dev_raw[mcquery_dev_raw[idx].originNumber * (2 * k) + j] = 				
					thrust::make_pair( (vtx_dev_raw[mcquery_dev_raw[idx].originNumber].r() - vtx_dev_raw[mcdata_dev_raw[initneig_dev_raw[idx] + j].originNumber].r()).length2(),
						(size_t)mcdata_dev_raw[initneig_dev_raw[idx] + j].originNumber );		
			}
		});


		thrust::for_each(
			thrust::device,
			thrust::make_counting_iterator<unsigned int>(0),
			thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
			[initdist_dev_raw, dstKeys1_dev_raw, k]
			__device__(unsigned int idx)
		{
			newSortcuda(initdist_dev_raw + (idx * 2 * k), initdist_dev_raw + (idx + 1) * 2 * k, dstKeys1_dev_raw + (idx * 2 * k));
		});
	}
	else
	{		
		thrust::pair<double, size_t>* dist_dev_raw = thrust::raw_pointer_cast(&dist_dev[0]);

		thrust::for_each(
			thrust::device,
			thrust::make_counting_iterator<unsigned int>(0),
			thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
			[dist_dev_raw, mcquery_dev_raw, mcdata_dev_raw, vtx_dev_raw, initneig_dev_raw, k, nvtx]
			__device__(unsigned int idx)
		{
				for (size_t j = 0; j < 2 * k; ++j)
				{
					dist_dev_raw[mcquery_dev_raw[idx].originNumber * 2 * k + j] =
					{ (vtx_dev_raw[mcquery_dev_raw[idx].originNumber].r() - vtx_dev_raw[mcdata_dev_raw[initneig_dev_raw[idx] + j].originNumber].r()).length2(),
					  (size_t)(mcdata_dev_raw[initneig_dev_raw[idx] + j].originNumber) };
				}

		});
				
		size_t* loc_dev_raw = thrust::raw_pointer_cast(&loc_dev[0]);
		size_t* counter_dev_raw = thrust::raw_pointer_cast(&counter_dev[0]);
		size_t* offset_dev_raw = thrust::raw_pointer_cast(&offset_dev[0]);
		size_t* counterScan_dev_raw = thrust::raw_pointer_cast(&counterScan_dev[0]);
		thrust::pair<double, size_t>* updateNN_dev_raw = thrust::raw_pointer_cast(&updateNN_dev[0]);

		thrust::for_each(
			thrust::device,
			thrust::make_counting_iterator<unsigned int>(0),
			thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
			[initdist_dev_raw, dstKeys_dev_raw, mcquery_dev_raw, dist_dev_raw, loc_dev_raw, counter_dev_raw, offset_dev_raw, counterScan_dev_raw, updateNN_dev_raw, k]
			__device__(unsigned int idx)
		{
			newMergecuda(k,
				initdist_dev_raw + (mcquery_dev_raw[idx].originNumber) * 2 * k,
				initdist_dev_raw + (mcquery_dev_raw[idx].originNumber) * 2 * k + k,
				dist_dev_raw + (mcquery_dev_raw[idx].originNumber) * 2 * k,
				loc_dev_raw + idx * 2 * k,
				counter_dev_raw + idx * 2 * k,
				offset_dev_raw + idx * 2 * k,
				counterScan_dev_raw + idx * 2 * k,
				updateNN_dev_raw + idx * k);
		});



		thrust::for_each(
			thrust::device,
			thrust::make_counting_iterator<unsigned int>(0),
			thrust::make_counting_iterator<unsigned int>((unsigned int)nvtx),
			[initdist_dev_raw, dstKeys_dev_raw, k]
			__device__(unsigned int idx)
		{
			newSortcuda(initdist_dev_raw + (idx * 2 * k), initdist_dev_raw + (idx) * 2 * k + k, dstKeys_dev_raw + (idx * k));
		});


	}
	//time[7] = omp_get_wtime();
}


void kNNcuda(const std::vector<Vortex2D>& vtx, 
	const size_t k, 
	std::vector<std::pair<double, size_t>>& initdist, 
	VectorsForKnn*& vecForKnn)
{
	//double t1 = omp_get_wtime();
	size_t nvtx = vtx.size();
	
	if (vecForKnn == nullptr)
		vecForKnn = new VectorsForKnn;

	vecForKnn->ResizeAll((int)k, (int)nvtx);
	
	Vortex2D* devVtxPtrV2D = thrust::raw_pointer_cast(&vecForKnn->vtx_dev[0]);
	cudaMemcpy(devVtxPtrV2D, vtx.data(), nvtx * sizeof(Vortex2D), cudaMemcpyHostToDevice);
	
	auto gab = GetGabaritcuda(devVtxPtrV2D, (int)nvtx);
	double scale = std::max(gab.second[0] - gab.first[0], gab.second[1] - gab.first[1]);
	//Сдвигаем все точки
	double time[4][8];
	for (size_t sdvig = 0; sdvig <= 4; ++sdvig) //4
	{

		//tStart[sdvig] = omp_get_wtime();
		mainCyclecuda((int)k, vtx, devVtxPtrV2D, gab, scale, (int)sdvig,
			vecForKnn->mcdatadev, vecForKnn->mcquerydev,
			vecForKnn->iqdev, vecForKnn->initneigdev, time[sdvig],
			vecForKnn->unsorted_mortonCodeBoth_dev, vecForKnn->mortonCodeBoth_dev,
			vecForKnn->unsorted_indicesBoth_dev, vecForKnn->indicesBoth_dev,
			vecForKnn->mortonCodeIndexBoth_dev,

			vecForKnn->prfx,
			vecForKnn->posTrue,
			vecForKnn->dstKeys1_dev,
			vecForKnn->dstKeys_dev,
			vecForKnn->initdist_dev,
			vecForKnn->dist_dev,
			vecForKnn->loc_dev,
			vecForKnn->counter_dev,
			vecForKnn->offset_dev,
			vecForKnn->counterScan_dev,
			vecForKnn->updateNN_dev
			);
		

		/*
		std::cout << "time_knn:" \
			<< " " << time[sdvig][1] - time[sdvig][0] \
			<< " " << time[sdvig][2] - time[sdvig][1] \
			<< " " << time[sdvig][3] - time[sdvig][2] \
			<< " " << time[sdvig][4] - time[sdvig][3] \
			<< " " << time[sdvig][5] - time[sdvig][4] \
			<< " " << time[sdvig][6] - time[sdvig][5] \
			<< " " << time[sdvig][7] - time[sdvig][6] << std::endl;
		*/
	}//for sdvig

	thrust::host_vector<thrust::pair<double, size_t>> initdist_host(vecForKnn->initdist_dev.begin(), vecForKnn->initdist_dev.begin() + vtx.size() * 2 * k);

#pragma omp parallel for
	for (int i = 0; i < vtx.size(); ++i)
	{
		for (int j = 0; j < k; ++j)
		{
			initdist[i * k + j].first = initdist_host[i * 2 * k + j].first;
			initdist[i * k + j].second = initdist_host[i * 2 * k + j].second;
		}
	}
}

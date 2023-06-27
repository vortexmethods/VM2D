
#include "knnCPU.h"

#include <omp.h>
#include <algorithm>
#include <iostream>
#include <bitset>
#include <tuple>


const int codeLength = 15;
const int twoPowCodeLength = (1 << codeLength);



//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
unsigned int ExpandBits(unsigned int v)
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
unsigned int Morton2D(const Point2D& r)
{
	const Point2D& rscale = twoPowCodeLength * r;
	const unsigned int& xx = ExpandBits((unsigned int)(rscale[0]));
	const unsigned int& yy = ExpandBits((unsigned int)(rscale[1]));
	return yy | (xx << 1);
}


/// Сортировка массива из мортоновский кодов 
void RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s)
{
	if (n == 0)
		return;
	// Количество задействованных потоков
	unsigned char threads = omp_get_max_threads();
	//std::cout << "threads = " << (int)threads << std::endl;

#pragma omp parallel num_threads(threads)
	{
		TParticleCode* source = m;
		TParticleCode* dest = m_temp;
		unsigned int l = omp_get_thread_num();
		unsigned int div = n / omp_get_num_threads();
		unsigned int mod = n % omp_get_num_threads();
		unsigned int left_index = l < mod ? (div + (mod == 0 ? 0 : 1)) * l : n - (omp_get_num_threads() - l) * div;
		unsigned int right_index = left_index + div - (mod > l ? 0 : 1);

		for (unsigned int digit = 0; digit < sizeof(m->key); ++digit)
		{
			unsigned int s_sum[256] = { 0 };
			unsigned int s0[256] = { 0 };
			unsigned char* b1 = (unsigned char*)&source[right_index].key;
			unsigned char* b2 = (unsigned char*)&source[left_index].key;
			while (b1 >= b2)
			{
				++s0[*(b1 + digit)];
				b1 -= sizeof(TParticleCode);
			}
			for (unsigned int i = 0; i < 256; i++)
			{
				s[i + 256 * l] = s0[i];
			}

#pragma omp barrier
			for (unsigned int j = 0; j < threads; j++)
			{
				for (unsigned int i = 0; i < 256; i++)
				{
					s_sum[i] += s[i + 256 * j];
					if (j < l)
					{
						s0[i] += s[i + 256 * j];
					}
				}
			}

			for (unsigned int i = 1; i < 256; ++i)
			{
				s_sum[i] += s_sum[i - 1];
				s0[i] += s_sum[i - 1];
			}
			unsigned char* b = (unsigned char*)&source[right_index].key + digit;
			TParticleCode* v1 = &source[right_index];
			TParticleCode* v2 = &source[left_index];
			while (v1 >= v2)
			{
				dest[--s0[*b]] = *v1--;
				b -= sizeof(TParticleCode);
			}
#pragma omp barrier
			std::swap(source, dest);
		}
	}

	// Если ключ структуры однобайтовый, просто копируем в исходный массив
	if (sizeof(m->key) == 1)
	{
		memcpy(m, m_temp, n * sizeof(TParticleCode));
	}
}


inline size_t BinSearch(const std::vector<std::pair<double, size_t>>& currentNN, double x, int low, int high)
{
	int mid = -1;

	if (x > currentNN[high].first)
		return high + 1;

	while (low <= high) {
		mid = (low + high) / 2;

		if (currentNN[mid].first == x)
			return mid + 1; //выход из цикла

		if (currentNN[mid].first < x)
			low = mid + 1;
		else
			high = mid - 1;
	}
	return mid;
}

void newSort(std::vector<std::pair<double, size_t>>& mass, std::vector<std::pair<double, size_t>>& dstKeys) {
	size_t k = mass.size();
	size_t cnt = 0;

	for (size_t i = 0; i < k; i++) {
		double elem = mass[i].first;

		cnt = 0;
		for (size_t j = 0; j < k; j++) {
			cnt += (mass[j].first < elem);
			//std::cout << cnt;
		}
		//if (!AssumeDistinct) 
		for (size_t j = 0; j < i; j++)
			cnt += (mass[j].first == elem);
		dstKeys[cnt] = mass[i];
	}
	mass.swap(dstKeys);
}

void newMerge(std::vector<std::pair<double, size_t>>& currentNN, const std::vector<std::pair<double, size_t>>& candidateNN,
	std::vector<size_t>& loc,
	std::vector<size_t>& counter,
	std::vector<size_t>& offset,
	std::vector<size_t>& counterScan,
	std::vector<std::pair<double, size_t>>& updateNN
)
{
	const size_t k = candidateNN.size() / 2;  //количество ближайших соседей

	//std::cout << "loc: " << std::endl;		
	for (size_t j = 0; j < 2 * k; ++j)
	{
		loc[j] = BinSearch(currentNN, candidateNN[j].first, 0, (int)k - 1);
		//std::cout << loc[j] << " ";
	}
	//std::cout << std::endl;

	for (size_t j = 0; j < 2 * k; ++j)
		counter[j] = j & 1;

	//std::cout << "counter: " << std::endl;
	//for (size_t j = 0; j < 2 * k; ++j)
	//    std::cout << counter[j] << " ";
	//std::cout << std::endl;

	for (size_t j = 0; j < 2 * k; ++j)
	{
		if ((loc[j] > 0) && ((loc[j] == k) || (candidateNN[j].second == currentNN[loc[j] - 1].second)))
			offset[j] = k + 1;
		else
			offset[j] = counter[2 * loc[j]]++;
	}

	//std::cout << "counter: " << std::endl;
	//for (size_t j = 0; j < 2 * k; ++j)
	//    std::cout << counter[j] << " ";
	//std::cout << std::endl;

	//std::cout << "offset: " << std::endl;
	//for (size_t j = 0; j < 2 * k; ++j)
	//    std::cout << offset[j] << " ";
	//std::cout << std::endl;


	counterScan[0] = 0;
	for (size_t j = 1; j < 2 * k; ++j)
		counterScan[j] = counterScan[j - 1] + counter[j - 1];

	//std::cout << "counterScan: " << std::endl;
	//for (size_t j = 0; j < 2 * k; ++j)
	//    std::cout << counterScan[j] << " ";
	//std::cout << std::endl;



	size_t index;
	for (size_t j = 0; j < k; ++j)
	{
		index = counterScan[2 * j + 1];
		//std::cout << index << std::endl;
		if (index < k)
			updateNN[index] = currentNN[j];
	}
	//std::cout << "updateNN: " << std::endl;
	//for (size_t j = 0; j < k; ++j)
	//    std::cout << "(" << updateNN[j].first << ", " << updateNN[j].second << "), ";
	//std::cout << std::endl;

	for (size_t j = 0; j < 2 * k; ++j)
	{
		if (2 * loc[j] < 2 * k)
		{
			index = counterScan[2 * loc[j]] + offset[j];
			if (index < k)
				updateNN[index] = candidateNN[j];
		}
	}

	currentNN.swap(updateNN);
}


std::pair<bool, bool> calcCheck(const Vortex2D& vtxi, const Vortex2D& vtxk, double maxG, double cSP, double cRBP, double epsCol, int type, double& d2)
{
	bool flagExit = false;
	bool check = false;

	//линейное увеличение радиуса коллапса				    
	double mnog = std::max(1.0, /* 2.0 * */ (vtxi.r()[0] - cRBP) / cSP);
	double r2test = (epsCol * mnog) * (epsCol * mnog);

	if (type == 1)
		r2test *= 4.0; //Увеличение радиуса коллапса в 2 раза для коллапса вихрей разных знаков		

	d2 = (vtxi.r() - vtxk.r()).length2();

	if (d2 > 1.0 * r2test)
	{
		flagExit = true;
		return { true, false };
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

	return { false, check };
}



void WakekNN(const std::vector<Vortex2D>& vtx, const size_t k, std::vector<std::vector<std::pair<double, size_t>>>& initdist,
	double cSP, double cRBP, double maxG, double epsCol, int type) 
{
	
	/////////////// ТЕСТ /////////////// 
	/*
	size_t Tk = 5;
	std::vector<std::pair<double, size_t>> Tinitdist{ {0.3, 1}, {2.3, 2}, {2.5, 333}, {5.4, 4}, {8.5, 5} };
	std::vector<std::pair<double, size_t>> TdstKeys1(Tk);
	std::vector<std::pair<double, size_t>> TdstKeys(2*Tk);


	std::vector<std::pair<double, size_t>> Tdist{ 
		{1.9, 101}, {5.0, 102}, {2.0, 103}, {8.7, 104}, {1.2, 105}, {9.0, 106}, {3.7, 107}, {2.5, 333}, {8.4, 109}, {9.3, 110} };
	std::vector<size_t> Tloc(2 * Tk);
	std::vector<size_t> Tcounter(2 * Tk);
	std::vector<size_t> Toffset(2 * Tk);
	std::vector<size_t> TcounterScan(2 * Tk);
	std::vector<std::pair<double, size_t>> TupdateNN(Tk);

	newSort(Tinitdist, TdstKeys1);
	Tinitdist.resize(Tk);
	newSort(Tdist, TdstKeys);
	newMerge(Tinitdist, Tdist, Tloc, Tcounter, Toffset, TcounterScan, TupdateNN);
	//newSort(Tinitdist, TdstKeys);
	*/	
	
	
	
	double preTime = -omp_get_wtime();
	
	const size_t n = vtx.size();

	/// Поиск границ габаритного прямоугольника (левый нижний и правый верхний) 
	auto xx = std::minmax_element(vtx.begin(), vtx.end(), [](const Vortex2D& P1, const Vortex2D& P2) {return P1.r()[0] < P2.r()[0]; });
	auto yy = std::minmax_element(vtx.begin(), vtx.end(), [](const Vortex2D& P1, const Vortex2D& P2) {return P1.r()[1] < P2.r()[1]; });
	
	double scale = std::max(xx.second->r()[0] - xx.first->r()[0], yy.second->r()[1] - yy.first->r()[1]);

	//создание единого массива data и query 
	std::vector<TParticleCode> mortonCodeBoth(2 * vtx.size());

	//сортировка единого массива (data и query) по мортоновским кодам
	std::vector<unsigned int> s(256 * omp_get_max_threads());
	std::vector<TParticleCode> mortonCodeBoth_temp(2 * vtx.size());

	//для каждого запроса запоминаем индекс iq -- индекс ближайшего левого элемента данных
	//отсортиванный массив mortonCodeBoth разделяем на mcquery (только запросы) и mcdata (только данные)
	std::vector<TParticleCode> mcdata, mcquery;
	std::vector<unsigned int> iq;

	//Сдвигаем все точки
	double t1 = omp_get_wtime();
	double tStart[5], tFinish[5];
	double tm[5][7];

	const int nSdvig = 3;

	preTime += omp_get_wtime();

	//std::cout << "preTime = " << preTime << std::endl;

	//14-05-2024
	for (size_t sdvig = 0; sdvig < nSdvig /*5*/; ++sdvig) 
	{
		tStart[sdvig] = omp_get_wtime();	
		
		const Point2D& lowLeft = { xx.first->r()[0], yy.first->r()[1] };
		double* time = tm[sdvig];		
		
		//масштабирование координат вихрей в [0;0.75)^2, поиск их мортоновских кодов 

		time[0] = omp_get_wtime();
#pragma omp parallel for
		for (int i = 0; i < vtx.size(); ++i)
		{
			//создание массива data 
			mortonCodeBoth[i].key = Morton2D((vtx[i].r() - lowLeft) * (0.75 / scale) + sdvig * Point2D{ 0.05, 0.05 });
			mortonCodeBoth[i].key <<= 1;

			mortonCodeBoth[i].originNumber = i;

			//создание массива query
			mortonCodeBoth[vtx.size() + i].key = mortonCodeBoth[i].key;

			mortonCodeBoth[vtx.size() + i].key |= 1;
			mortonCodeBoth[vtx.size() + i].originNumber = i;
		}

		time[1] = omp_get_wtime();
	
		//сортировка единого массива (data и query) по мортоновским кодам			
		RSort_Parallel(mortonCodeBoth.data(), mortonCodeBoth_temp.data(), (int)mortonCodeBoth.size(), s.data());

		time[2] = omp_get_wtime();

		//для каждого запроса запоминаем индекс iq -- индекс ближайшего левого элемента данных
		//отсортиванный массив mortonCodeBoth разделяем на mcquery (только запросы) и mcdata (только данные)
		mcdata.resize(0);
		mcquery.resize(0);
		mcdata.reserve(vtx.size());
		mcquery.reserve(vtx.size());

		time[3] = omp_get_wtime();

		iq.resize(0);
		iq.reserve(vtx.size());
		for (size_t i = 0, j = -1; i < mortonCodeBoth.size(); ++i)
		{
			if (mortonCodeBoth[i].key & 1)
			{
				mcquery.push_back(mortonCodeBoth[i]);
				iq.push_back((unsigned)(j + 1) == 0 ? 0 : (unsigned int)j);
			}
		else
			{
				++j;
				mcdata.push_back(mortonCodeBoth[i]);
			}
		}

		time[4] = omp_get_wtime();
		
		std::vector<std::pair<double, size_t>> dist(2 * k, { -1.0, 0 });

		std::vector<size_t> loc(2 * k);
		std::vector<size_t> counter(2 * k, 0);
		std::vector<size_t> offset(2 * k, 0);
		std::vector<size_t> counterScan(2 * k, 0);
		std::vector<std::pair<double, size_t>> updateNN(k, { 0.0, 0 });

		std::vector<std::pair<double, size_t>> dstKeys1(2 * k, { 0.0, 0 });
		std::vector<std::pair<double, size_t>> dstKeys(k, { 0.0, 0 });

		time[5] = omp_get_wtime();

//*
#pragma omp parallel for firstprivate(dist, loc, counter, offset, counterScan, updateNN, dstKeys, dstKeys1)
		for (int i = 0; i < initdist.size(); ++i)
		{
			const Vortex2D& vtxi = vtx[mcquery[i].originNumber];

			int cntr = 0;
			int search = iq[i];
			std::vector<std::pair<double, size_t>>& fillPosition = ((sdvig == 0) ? initdist[mcquery[i].originNumber] : dist);

			while ((cntr < k) && (search >= 0))
			{
				const Vortex2D& vtxk = vtx[mcdata[search].originNumber];
				bool flagExit, check;
				double d2;
				std::tie(flagExit, check) = calcCheck(vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2);

				if (flagExit)
					break;

				if ((mcdata[search].originNumber > mcquery[i].originNumber) && check)
					fillPosition[cntr++] = { d2, mcdata[search].originNumber };

				--search;
			}

			for (int w = cntr; w < k; ++w)
				fillPosition[cntr++] = { 100000000.0, 0 };

			search = iq[i] + 1;

			while ((cntr < 2 * k) && (search < mcdata.size()))
			{
				const Vortex2D& vtxk = vtx[mcdata[search].originNumber];
				bool flagExit, check;
				double d2;
				std::tie(flagExit, check) = calcCheck(vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2);

				if (flagExit)
					break;

				if ((mcdata[search].originNumber > mcquery[i].originNumber) && check)
					fillPosition[cntr++] = { d2,  mcdata[search].originNumber };

				++search;
			}
			for (int w = cntr; w < 2 * k; ++w)
				fillPosition[cntr++] = { 100000000.0, 0 };






			if (sdvig == 0)
			{
				newSort(initdist[mcquery[i].originNumber], dstKeys1);
				initdist[mcquery[i].originNumber].resize(k);
			}
			else
			{
				newSort(dist, dstKeys1);
				newMerge(initdist[mcquery[i].originNumber], dist, loc, counter, offset, counterScan, updateNN);
				//newSort(initdist[mcquery[i].originNumber], dstKeys);
			}
		}
//*/

		time[6] = omp_get_wtime();
		
		tFinish[sdvig] = omp_get_wtime();

		//std::cout << "time_knn:" 
		//	<< " " << tm[sdvig][1] - tm[sdvig][0] 
		//	<< " " << tm[sdvig][2] - tm[sdvig][1] 
		//	<< " " << tm[sdvig][3] - tm[sdvig][2] 
		//	<< " " << tm[sdvig][4] - tm[sdvig][3] //0.0002688
		//	<< " " << tm[sdvig][5] - tm[sdvig][4] 
		//	<< " " << tm[sdvig][6] - tm[sdvig][5] //0.0046687
		//  << std::endl;
	}//for sdvig
}
#ifndef WRAPPER_H_
#define WRAPPER_H_

#include "cuKernels.cuh"

namespace BHcu
{
	double memoryAllocate(CUDApointers& ptrs, int nnodes, int nbodies, int nbodiesOld, int blocks, int order);
	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="vtxl"></param> ”казатель на вихри на GPU
	/// <param name="vell"></param> ”казатель на массив вычисл€емых скоростей на GPU
	/// <param name="epsastl"></param> ”казатель на массив eps* на GPU
	/// <param name="ptrs"></param> ’ранилище указателей, относ€щихс€ к дереву
	/// <param name="nbodies"></param> число вихрей
	/// <param name="timing"></param> указатель на 7-элементный массив дл€ засечек времени дл€ быстрого метода
	/// <param name="eps"></param> eps
	/// <param name="theta"></param> theta
	/// <param name="nbodiesOld"></param> число вихрей, под которые была выделена пам€ть на прошлом шаге
	/// <param name="nbodiesUp"></param> число вихрей (оценка сверху), под которые требуетс€ пам€ть на текущем шаге
	/// <param name="order"></param> order
	/// <param name="nAfls"></param> число профилей
	/// <param name="nVtxs"></param> указатель на массив из числа вихрей, сгенерированных на профил€ъ
	/// <param name="ptrVtxs"></param> указатели на них
	/// <returns></returns>
	
	double wrapperInfluence(const realVortex* vtxl, realPoint* vell, real* epsastl, CUDApointers& ptrs, int nbodies, 
		double* timing, real eps, real theta, size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs);

	/// npoints - число точек наблюдени€
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild,		
		int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs);


	double wrapperInfluenceToRHS(
		const realVortex* dev_ptr_vt,  //вихри в следе
		const double* dev_ptr_pt,      //начала и концы панелей
		double* dev_ptr_rhs,           //куда сохранить результат 
		double* dev_ptr_rhslin,        //куда сохранить результат 

		CUDApointers& ptrs,     //указатели на делево вихрей
		bool rebuild,           //признак перестроени€ делева вихрей

		int nvt,               //число вихрей в следе
		int nTotPan,           //общее число панелей на всех профил€х
		double* timingsToRHS,  //засечки времени
		//double eps,            //eps
		double theta,          //theta
		size_t& nbodiesOld, int nbodiesUp, int order, int scheme);

	double wrapperDiffusiveVelo(const realVortex* vtxl,
		real* i1l,
		realPoint* i2l,
		real* epsastl,
		CUDApointers& ptrs,
		bool rebuild,
		int nbodies,
		double* timing,
		real eps, real theta,
		size_t& nbodiesOld,
		int nbodiesUp,
		int order,
		size_t nAfls,
		size_t* nVtxs,
		double** ptrVtxs);

}

#endif
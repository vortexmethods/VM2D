#ifndef WRAPPER_H_
#define WRAPPER_H_

#include "cuKernels.cuh"

namespace BHcu
{

//*
	class TMortonCodesCalculator
	{
	protected:
		int nPoints;
		const void* ptr;
		int sizeOfArrayElement;
		int offsetOfPointInElement;

		int blocks;
		bool sort;	

		CUDApointers cudaPtrs;

	public:
		const int* getSortedIndices() const { return cudaPtrs.MmortonCodesIdxl; }

		float timeOfWork;

		TMortonCodesCalculator(int nPoints_, const void* ptr_, int sizeOfArrayElement_, int offsetOfPointInElement_, bool sort_ = true)
			: nPoints(nPoints_), ptr(ptr_), sizeOfArrayElement(sizeOfArrayElement_), offsetOfPointInElement(offsetOfPointInElement_), sort(sort_)
		{
			cudaPtrs.MmortonCodesKeyUnsortl = (int*)cudaNew(nPoints, sizeof(int));
			cudaPtrs.MmortonCodesIdxUnsortl = (int*)cudaNew(nPoints, sizeof(int));

			if (sort)
			{
				cudaPtrs.MmortonCodesKeyl = (int*)cudaNew(nPoints, sizeof(int));
				cudaPtrs.MmortonCodesIdxl = (int*)cudaNew(nPoints, sizeof(int));
			}			
				
			CudaSelect(0);
			setBlocks(blocks);

			cudaPtrs.maxrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
			cudaPtrs.minrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));			

			float timeBoundingBox = McuBoundingBoxKernelFree(nullptr, cudaPtrs.maxrl, cudaPtrs.minrl, nPoints, ptr, sizeOfArrayElement, offsetOfPointInElement);

			float timeCodes = McuMortonCodesKernelFree(cudaPtrs.maxrl, cudaPtrs.minrl, cudaPtrs.MmortonCodesKeyUnsortl, cudaPtrs.MmortonCodesIdxUnsortl, cudaPtrs.MmortonCodesKeyl, cudaPtrs.MmortonCodesIdxl, nullptr, nPoints_,
				ptr, sizeOfArrayElement, offsetOfPointInElement, sort);

			timeOfWork = timeBoundingBox + timeCodes;
		};

		~TMortonCodesCalculator()
		{
			cudaDelete(cudaPtrs.MmortonCodesKeyUnsortl, 1);
			cudaDelete(cudaPtrs.MmortonCodesIdxUnsortl, 3);
			
			if (sort)
			{
				cudaDelete(cudaPtrs.MmortonCodesKeyl, 0);
				cudaDelete(cudaPtrs.MmortonCodesIdxl, 2);
			}

			cudaDelete(cudaPtrs.maxrl, 4);
			cudaDelete(cudaPtrs.minrl, 5);
		}
	};
//*/





	double memoryAllocate(CUDApointers& ptrs, int nnodes, int nbodies, int nbodiesOld, int blocks, int order);
	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="vtxl"></param> Указатель на вихри на GPU
	/// <param name="vell"></param> Указатель на массив вычисляемых скоростей на GPU
	/// <param name="epsastl"></param> Указатель на массив eps* на GPU
	/// <param name="ptrs"></param> Хранилище указателей, относящихся к дереву
	/// <param name="nbodies"></param> число вихрей
	/// <param name="timing"></param> указатель на 7-элементный массив для засечек времени для быстрого метода
	/// <param name="eps"></param> eps
	/// <param name="theta"></param> theta
	/// <param name="nbodiesOld"></param> число вихрей, под которые была выделена память на прошлом шаге
	/// <param name="nbodiesUp"></param> число вихрей (оценка сверху), под которые требуется память на текущем шаге
	/// <param name="order"></param> order
	/// <param name="nAfls"></param> число профилей
	/// <param name="nVtxs"></param> указатель на массив из числа вихрей, сгенерированных на профиляъ
	/// <param name="ptrVtxs"></param> указатели на них
	/// <returns></returns>
	/// npoints - число точек наблюдения
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, int sizeOfArrayElement, int offsetOfPointInElement,
		realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild,		
		int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, 
		int nbodiesUp, 
		int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs, bool calcVelo, bool calcRadius);


	double wrapperInfluenceToRHS(
		const realVortex* dev_ptr_vt,  //вихри в следе
		const double* dev_ptr_pt,      //начала и концы панелей
		double* dev_ptr_rhs,           //куда сохранить результат 
		double* dev_ptr_rhslin,        //куда сохранить результат 

		CUDApointers& ptrs,     //указатели на делево вихрей
		bool rebuild,           //признак перестроения делева вихрей

		int nvt,               //число вихрей в следе
		int nTotPan,           //общее число панелей на всех профилях
		double* timingsToRHS,  //засечки времени
		//double eps,            //eps
		double theta,          //theta
		size_t& nbodiesOld, int nbodiesUp, int order, int scheme);


	class wrapperMatrixToVector
	{
		const double* dev_ptr_pt;      //начала и концы панелей
		double* dev_ptr_rhs;          //куда сохранить результат (для T0 и верхней половины T1)
		double* dev_ptr_rhslin;        //куда сохранить результат (для нижней половины T1)

		CUDApointers& ptrs;     //указатели на дерево вихрей
		bool rebuild;         //признак перестроения дерева вихрей

		int nTotPan;           //общее число панелей на всех профилях
		double* timingsMatr;  //засечки времени
		double theta;          //theta
		int order; 
		int scheme;

		double* panelPoints;
		Point2D* controlPoints;
		realPoint* El;
		TMortonCodesCalculator* mCodesPtr;

	public:
		wrapperMatrixToVector(
			const double* dev_ptr_pt_,      //начала и концы панелей
			double* dev_ptr_rhs_,           //куда сохранить результат (для T0 и верхней половины T1)
			double* dev_ptr_rhslin_,        //куда сохранить результат (для нижней половины T1)

			CUDApointers& ptrs_,     //указатели на дерево вихрей
			bool rebuild_,           //признак перестроения дерева вихрей

			int nTotPan_,           //общее число панелей на всех профилях
			double* timingsMatr_,  //засечки времени
			double theta_,          //theta
			int order_, int scheme_);

		~wrapperMatrixToVector();

		double calculate(
			double* dev_ptr_freeVortexSheet,
			double* dev_ptr_freeVortexSheetLin);
	};
	 



	double wrapperInfluenceFromPanelsToPoints(
		double* dev_ptr_r,  //начала и концы панелей
		double* dev_ptr_freeVortexSheet, double* dev_ptr_freeVortexSheetLin,
		double* dev_ptr_attachedVortexSheet, double* dev_ptr_attachedVortexSheetLin,
		double* dev_ptr_attachedSourceSheet, double* dev_ptr_attachedSourceSheetLin,
		Vortex2D* dev_ptr_pt, //вихри в следе
		double*& dev_ptr_vel,  //куда сохранить результат 
		CUDApointers& ptrsi,  //указатель на дерево на i-м профиле
		bool rebuild,         //признак перестроения дерева вихрей
		int npt,			  //число вихрей в следе
		int npnli,			  //общее число панелей на профиле
		double* timings,      //засечки времени
		double eps,           //eps
		double theta,         //theta
		int order,            //order
		int schemeType
	);



	double wrapperDiffusiveVeloI1I2(const realVortex* vtxl,
		real* i1l,
		realPoint* i2l,
		real* epsastl,
		CUDApointers& ptrs,
		bool rebuild,
		int nbodies,
		double* timing,
		real minRd,
		size_t& nbodiesOld,
		int nbodiesUp,
		int order);

	double wrapperDiffusiveVeloI0I3(const realVortex* vtxl, 
		float* i0l, 
		floatPoint* i3l,
		real* epsastl, 
		real* ptr_r, 
		CUDApointers& ptrs, 
		bool rebuild, 
		int nbodies, 
		size_t nPan,
		real* visstr,
		double* timing, 
		real* meanEps, 
		real minRd, 
		size_t& nbodiesOld, 
		int nbodiesUp, 
		int order
		);





}

#endif
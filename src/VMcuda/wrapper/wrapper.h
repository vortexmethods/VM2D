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
	/// <param name="vtxl"></param> ��������� �� ����� �� GPU
	/// <param name="vell"></param> ��������� �� ������ ����������� ��������� �� GPU
	/// <param name="epsastl"></param> ��������� �� ������ eps* �� GPU
	/// <param name="ptrs"></param> ��������� ����������, ����������� � ������
	/// <param name="nbodies"></param> ����� ������
	/// <param name="timing"></param> ��������� �� 7-���������� ������ ��� ������� ������� ��� �������� ������
	/// <param name="eps"></param> eps
	/// <param name="theta"></param> theta
	/// <param name="nbodiesOld"></param> ����� ������, ��� ������� ���� �������� ������ �� ������� ����
	/// <param name="nbodiesUp"></param> ����� ������ (������ ������), ��� ������� ��������� ������ �� ������� ����
	/// <param name="order"></param> order
	/// <param name="nAfls"></param> ����� ��������
	/// <param name="nVtxs"></param> ��������� �� ������ �� ����� ������, ��������������� �� ��������
	/// <param name="ptrVtxs"></param> ��������� �� ���
	/// <returns></returns>
	/// npoints - ����� ����� ����������
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, int sizeOfArrayElement, int offsetOfPointInElement,
		realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild,		
		int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs, bool calcVelo, bool calcRadius);


	double wrapperInfluenceToRHS(
		const realVortex* dev_ptr_vt,  //����� � �����
		const double* dev_ptr_pt,      //������ � ����� �������
		double* dev_ptr_rhs,           //���� ��������� ��������� 
		double* dev_ptr_rhslin,        //���� ��������� ��������� 

		CUDApointers& ptrs,     //��������� �� ������ ������
		bool rebuild,           //������� ������������ ������ ������

		int nvt,               //����� ������ � �����
		int nTotPan,           //����� ����� ������� �� ���� ��������
		double* timingsToRHS,  //������� �������
		//double eps,            //eps
		double theta,          //theta
		size_t& nbodiesOld, int nbodiesUp, int order, int scheme);

	double wrapperInfluenceFromPanelsToPoints(
		double* dev_ptr_r,  //������ � ����� �������
		double* dev_ptr_freeVortexSheet, double* dev_ptr_freeVortexSheetLin,
		double* dev_ptr_attachedVortexSheet, double* dev_ptr_attachedVortexSheetLin,
		double* dev_ptr_attachedSourceSheet, double* dev_ptr_attachedSourceSheetLin,
		Vortex2D* dev_ptr_pt, //����� � �����
		double*& dev_ptr_vel,  //���� ��������� ��������� 
		CUDApointers& ptrsi,  //��������� �� ������ �� i-� �������
		bool rebuild,         //������� ������������ ������ ������
		int npt,			  //����� ������ � �����
		int npnli,			  //����� ����� ������� �� �������
		double* timings,      //������� �������
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
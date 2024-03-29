#ifndef WRAPPER_H_
#define WRAPPER_H_

#include "cuKernels.cuh"

namespace BHcu
{
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
	
	double wrapperInfluence(const realVortex* vtxl, realPoint* vell, real* epsastl, CUDApointers& ptrs, int nbodies, 
		double* timing, real eps, real theta, size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs);

	/// npoints - ����� ����� ����������
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild,		
		int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs);


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
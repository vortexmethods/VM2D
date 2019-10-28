/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidRotateMon.cpp                                    |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM is distributed in the hope that it will be useful, but WITHOUT           |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса MechanicsRigidRotateMon
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/

#include "mpi.h"

#include "Mechanics2DRigidRotateMon.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidRotateMon::GetHydroDynamForce(timePeriod& time)
{
	time.first = omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };		//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	Point2D hDFQ = { 0.0, 0.0 };		//гидродинамические силы, обусловленные Q_k
	double hDMGam = 0.0;				//гидродинамический момент, обусловленный Gamma_k
	double hDMdelta = 0.0;				//гидродинамический момент, обусловленный Gamma_k
	double hDMQ = 0.0;					//гидродинамический момент, обусловленный Gamma_k

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;
		Point2D velK = { -Wcm * rK[1], Wcm * rK[0] };

		double gAtt = Wcm * (rK ^ afl.tau[i]);
		double gAttOld = WcmOld * (rK ^ afl.tau[i]);
		double deltaGAtt = gAtt - gAttOld;

		double qAtt = Wcm * (rK ^ afl.nrm[i]);
		
		/*1*/
		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * (rK * rK);
		
		/*2*/
		hDFGam += 0.5 * afl.getV(i).kcross() * gAtt * afl.len[i];
		hDMGam += 0.5 * rK ^ afl.getV(i).kcross() * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.5 * afl.getV(i) * qAtt * afl.len[i];
		hDMQ -= 0.5 * rK ^ afl.getV(i) * qAtt * afl.len[i];
	}

	hydroDynamForce = hDFGam + hDFdelta * (1.0 / dt) + hDFQ;
	hydroDynamMoment = hDMGam + hDMdelta / dt + hDMQ;

	time.second = omp_get_wtime();
}// CalcHydroDynamForce(...)

// Вычисление скорости центра масс
Point2D MechanicsRigidRotateMon::VeloOfAirfoilRcm(double currTime)
{
	return{ 0.0, 0.0 };
	//return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidRotateMon::PositionOfAirfoilRcm(double currTime)
{
	return afl.rcm;
}//PositionOfAirfoilRcm(...)

// Вычисление скоростей начал панелей
void MechanicsRigidRotateMon::VeloOfAirfoilPanels(double currTime)
{
	std::vector<Point2D> vel(afl.getNumberOfPanels(), { 0.0, 0.0 });
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		vel[i][0] = afl.getR(i).kcross()[0] * Wcm;
		vel[i][1] = afl.getR(i).kcross()[1] * Wcm;
	}
	afl.setV(vel);

}//VeloOfAirfoilPanels(...)

//Перемещение профиля в соответствии с законом
void MechanicsRigidRotateMon::Move()
{
	PhiOld = Phi;

	double dphi = Wcm * W.getPassport().timeDiscretizationProperties.dt;
	MPI_Bcast(&dphi, 1, MPI_DOUBLE, 0, W.getParallel().commWork);

	Phi += dphi;
	afl.Rotate(dphi);
}//Move()

//Чтение параметров конкретной механической системы
void MechanicsRigidRotateMon::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("J", J);
	mechParamsParser->get("k", k);
	mechParamsParser->get("tRotateAccel", tRotateAccel);
	mechParamsParser->get("tMomentAccel", tMomentAccel);
	mechParamsParser->get("Mz", Mz);

	W.getInfo('i') << "J = " << J << std::endl;
	W.getInfo('i') << "k = " << k << std::endl;
	W.getInfo('i') << "tRotateAccel = " << tRotateAccel << std::endl;
	W.getInfo('i') << "tMomentAccel = " << tMomentAccel << std::endl;
	W.getInfo('i') << "Mz = " << Mz << std::endl;

//	std::vector<double> sh;
//	mechParamsParser->get("sh", sh);
//	k = 0.0 * m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.V0().length2();
}

//Заполнение строк в матрице, соответствующих механической системе.А также пересечений строк и столбцов, соответствующих мех.системе
void MechanicsRigidRotateMon::FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross)
{
	if (W.getCurrentStep() == 0)
	{
		cross(0, 0) = -(/*3*/ J + /*4*/ W.getPassport().timeDiscretizationProperties.dt * b);

		Point2D riC;
#pragma warning (push)
#pragma warning (disable: 4101)
		double crossOmp;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(row, cross) private(riC, crossOmp)
		for (int loc_i = 0; loc_i < afl.getNumberOfPanels(); ++loc_i)
		{
			riC = 0.5 * (afl.getR(loc_i + 1) + afl.getR(loc_i)) - afl.rcm;
			crossOmp = /*2*/ 0.5 * riC * riC * (riC ^ afl.tau[loc_i]) * afl.len[loc_i];

			row(0, loc_i) = /*1*/ 0.5 * riC * riC * afl.len[loc_i];
#pragma omp atomic
			cross(0, 0) += crossOmp;
		}//pragma omp for
	}//if (W.getCurrentStep() == 0)
}//FillMechanicsRowsAndCols(...)


//Заполнение правых частей, соответствующих механической системе
void MechanicsRigidRotateMon::FillMechanicsRhs(std::vector<double>& rhs)
{
	rhs[0] = /*7*/-J * Wcm + /*8*/ W.getPassport().timeDiscretizationProperties.dt * k * Phi -\
				/*11*/GetMz(W.getCurrentStep() * W.getPassport().timeDiscretizationProperties.dt)* W.getPassport().timeDiscretizationProperties.dt;
	
	Point2D riC;

#pragma warning (push)
#pragma warning (disable: 4101)
	double rhsOmp;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(rhs) private(riC, rhsOmp)
	for (int loc_i = 0; loc_i < afl.getNumberOfPanels(); ++loc_i)
	{
		riC = 0.5 * (afl.getR(loc_i + 1) + afl.getR(loc_i)) - afl.rcm;

		rhsOmp = 0.0;
		rhsOmp += /*5*/ 0.5 * riC * riC * afl.gammaThrough[loc_i];
		rhsOmp += /*6*/ 0.5 * riC * riC * Wcm * (riC ^ afl.tau[loc_i]) * afl.len[loc_i];
		rhsOmp += /*9*/ 0.5 * afl.len[loc_i] * (Wcm * riC ^ afl.tau[loc_i]) * (-Wcm * riC) ^ riC * W.getPassport().timeDiscretizationProperties.dt;
		rhsOmp -= /*10*/ -0.5 * Wcm * riC * riC * afl.len[loc_i] * (riC ^ afl.nrm[loc_i]) * Wcm * W.getPassport().timeDiscretizationProperties.dt;

#pragma omp atomic
		rhs[0] += rhsOmp;
	}
}//FillMechanicsRhs(...)

//Заполнение правых частей (расщепленная схема) либо коэффициентов в матрице при скорости (монолитный подход), соответствующих присоединенным слоям
void MechanicsRigidRotateMon::FillAtt(Eigen::MatrixXd& col, Eigen::MatrixXd& rhs)
{
	W.getInfo('t') << "FillAtt: " << W.getParallel().myidWork << std::endl;

	if (W.getCurrentStep() == 0)
	{
		Point2D riC;


		std::vector<double> colMPI;
		colMPI.resize(col.rows(), 0.0);

		int id = W.getParallel().myidWork;
		
		/// \todo Почему в других мех.подсистемах нет обращения к инфраструктуре MPI?
		VMlib::parProp par = W.getParallel().SplitMPI(col.rows());

		std::vector<double> locColMPI;
		locColMPI.resize(par.myLen);


#pragma omp parallel for default(none) shared(locColMPI, par) private(riC)
		for (int locI = 0; locI < par.myLen; ++locI)
		{
			int i = locI + par.myDisp;

			riC = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;

			locColMPI[locI] = /*12 + 91*/ -0.5 * (riC ^ afl.tau[i]);

			for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
			{
				locColMPI[locI] += /*92*/ IDPI * (Alpha(riC - afl.getR(j + 1), riC - afl.getR(j)) * afl.tau[j]
					- Lambda(riC - afl.getR(j + 1), riC - afl.getR(j)) * afl.nrm[j]
					)* afl.tau[i] * (riC ^ afl.tau[j]);

				locColMPI[locI] += /*93*/ IDPI * (Alpha(riC - afl.getR(j + 1), riC - afl.getR(j)) * afl.nrm[j] +
					Lambda(riC - afl.getR(j + 1), riC - afl.getR(j)) * afl.tau[j]
					)* afl.tau[i] * (riC ^ afl.nrm[j]);
			}
		}//pragma omp for
		MPI_Gatherv(locColMPI.data(), par.myLen, MPI_DOUBLE, colMPI.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

		if (id == 0)
		for (size_t i = 0; i < (size_t)(col.rows()); ++i)
		{
			col(i, 0) = colMPI[i];
		}
	}//if (W.getCurrentStep() == 0)	
}//FillAtt(...)

//Извлечение из решения параметров механической системы
void MechanicsRigidRotateMon::SolutionToMechanicalSystem(Eigen::VectorXd& col)
{
	WcmOld = Wcm;

	Wcm = col(0);

	std::vector<Point2D> vel(afl.getNumberOfPanels(), { 0.0, 0.0 });
	for (size_t i = 0; i < afl.getNumberOfPanels() + 1; ++i)
		vel[i] = afl.getR(i).kcross() * Wcm;

	afl.setV(vel);
}//SolutionToMechanicalSystem(...)
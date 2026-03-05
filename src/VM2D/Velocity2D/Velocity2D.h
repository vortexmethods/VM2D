/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Velocity2D.h                                                     |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса Velocity
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef VELOCITY2D_H
#define VELOCITY2D_H

#include <memory>

#include "defs.h"
#include "cudaTreeInfo.h"

namespace VM2D
{

	class Airfoil;
	class Boundary;
	class WakeDataBase;
	class World2D;
	

	/*!
	\brief Структура, определяющая параметры виртуальных вихрей для отдельного профиля

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	struct VortexesParams
	{
		/// Вектор конвективных скоростей вихрей
		std::vector<Point2D> convVelo;

		/// Вектор диффузионных скоростей вихрей
		std::vector<Point2D> diffVelo;

		/// Вектор знаменателей (I0) диффузионных скоростей вихрей (обусловленных профилем)
		std::vector<double> I0;

		/// Вектор знаменателей (I1) диффузионных скоростей вихрей (обусловленных завихренностью)
		std::vector<double> I1;

		/// Вектор числителей (I2) диффузионных скоростей вихрей (обусловленных завихренностью)
		std::vector<Point2D> I2;

		/// Вектор числителей (I3) диффузионных скоростей вихрей (обусловленных профилем)
		std::vector<Point2D> I3;

		/// Вектор характерных радиусов вихревых доменов (eps*)
		std::vector<double> epsastWake;
	};

	/*!
	\brief Абстрактный класс, определяющий способ вычисления скоростей

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class Velocity
	{
	protected:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

	public:
		/// Струтура, определяющая параметры вихрей в следе
		VortexesParams wakeVortexesParams;

		/// Вектор струтур, определяющий параметры виртуальных вихрей для профилей
		std::vector<VortexesParams> virtualVortexesParams;


		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу 	
		Velocity(const World2D& W_) :
			W(W_)
		{
			virtualVortexesParams.resize(0);
		};

		/// \brief Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
		///
		/// \warning скорости приплюсовываются к уже имеющимся
		void CalcConvVelo();

		/// \brief Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек от следа
		///
		/// \param[in] pointsDb константная ссылка на базу данных пелены из вихрей, в которых надо сосчитать конвективные скорости
		/// \param[out] velo ссылка на вектор скоростей в требуемых точках
		/// \param[out] domainRadius ссылка на вектор радиусов вихревых доменов
		/// \param[in] calcVelo признак вычисления скоростей в точках
		/// \param[in] calcRadius признак вычисления радиусов доменов
		virtual void CalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius) = 0;		
#if defined (USE_CUDA)
		virtual void GPUCalcConvVeloToSetOfPointsFromWake(std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree, const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius) = 0;
		virtual void GPUCalcConvVelocityToSetOfPointsFromSheets(std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree, const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const = 0;
#endif

		/// Вычисление скоростей в точках wakeVP
		virtual void CalcVeloToWakeVP();

		/// \brief Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
		///
		/// \param[in] pointsDb константная ссылка на базу данных пелены из вихрей, в которых надо сосчитать диффузионные скорости
		/// \param[in] domainRadius константная ссылка на вектор радиусов вихревых доменов
		/// \param[in] vorticesDb константная ссылка на на базу данных пелены из вихрей,от которых надо сосчитать влияния на points
		/// \param[out] I1 ссылка на вектор величин I1 (знаменателей в диффузионных скоростях) в требуемых точках
		/// \param[out] I2 ссылка на вектор величин I2 (числителей в диффузионных скоростях) в требуемых точках
		void CalcDiffVeloI1I2ToSetOfPointsFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2);
		void CalcDiffVeloI1I2ToSetOfPointsFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2);
		
		void CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2);
		void CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2);
#if defined(USE_CUDA)
		void GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh = false);
		void GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh = false);
		void GPUDiffVeloFAST(const std::vector<double>& domainRadius, std::vector<double>& I1, std::vector<Point2D>& I2);
#endif

		/// \brief Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
		///
		/// Вызывает 4 раза функцию CalcDiffVeloToSetOfPoints
		///
		/// \warning скорости приплюсовываются к уже имеющимся
		void CalcDiffVeloI1I2();
		void CalcDiffVeloI0I3();


		/// \brief Контроль больших значений диффузионных скоростей
		/// 
		/// \param[in, out] diffVel ссылка на вектор диффузионных скоростей
		void LimitDiffVelo(std::vector<Point2D>& diffVel);

		/// \brief Вычисление диффузионных скоростей
		///
		/// Вызывается в CalcVortexVelo()
		void CalcDiffVelo();

		/// \brief Очистка старых массивов под хранение скоростей, выделение новой памяти и обнуление
		///
		/// Вызывается в CalcVortexVelo() на каждом шаге расчета перед непосредственно расчетом скоростей
		void ResizeAndZero();

		/// \brief Сохранение вязких напряжений
		///
		/// Вызывается в CalcDiffVelo()
		void SaveVisStress();


		/// \brief Генерация вектора влияния вихревого следа на профиль
		///
		/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
		/// 
		/// \param[in] afl константная ссылка на профиль, правая часть для которого вычисляется
		/// \param[out] wakeRhs ссылка на вектор влияния вихревого следа на ОДИН профиль
		void GetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const;
#if defined(USE_CUDA)
		void GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const;
		void GPUFASTGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const;
#endif
		void FillRhs(/*Eigen::VectorXd& rhs, Eigen::VectorXd& rhsSkos,*/ Eigen::VectorXd& rhsReord) const;


		/// Деструктор
		virtual ~Velocity() { };
	};

}//namespace VM2D

#endif

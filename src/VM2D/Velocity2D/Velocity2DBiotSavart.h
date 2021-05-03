/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Velocity2DBiotSavart.h                                           |
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
\brief Заголовочный файл с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef VELOCITYBIOTSAVART_H
#define VELOCITYBIOTSAVART_H

#include "Velocity2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий способ вычисления скоростей

	Способ вычисления скоростей
	- напрямую по закону Био --- Савара

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.10
	\date 17 мая 2021 г.
	*/
	class VelocityBiotSavart : public Velocity
	{
	public:
		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу 
		VelocityBiotSavart(const World2D& W_);

		/// Деструктор
		virtual ~VelocityBiotSavart();


		/// Вычисление скоростей в точках wakeVP
		void CalcVeloToWakeVP();

		/// \brief Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек от следа
		///
		/// \param[in] pointsDb константная ссылка на базу данных пелены из вихрей, в которых надо сосчитать конвективные скорости
		/// \param[out] velo ссылка на вектор скоростей в требуемых точках
		/// \param[out] domainRadius ссылка на вектор радиусов вихревых доменов
		/// \param[in] calcVelo признак вычисления скоростей в точках
		/// \param[in] calcRadius признак вычисления радиусов доменов
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		void CalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius);
#if defined (USE_CUDA)
		void GPUCalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius);
#endif	

		//реализация виртуальных функций
		virtual void CalcConvVelo() override;
		virtual void FillRhs(Eigen::VectorXd& rhs) const override;


		/// \brief Генерация вектора влияния вихревого следа на профиль
		///
		/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
		/// 
		/// \param[in] afl константная ссылка на профиль, правая часть для которого вычисляется
		/// \param[out] wakeRhs ссылка на вектор влияния вихревого следа на ОДИН профиль
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		void GetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const;
#if defined(USE_CUDA)
		void GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const;
#endif

		void CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2) override;
		void CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2) override;

	};

}//namespace VM2D

#endif
/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Airfoil2DRect.h                                                  |
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
\brief Заголовочный файл с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#ifndef AIRFOILCURV_H
#define AIRFOILCURV_H

#include "Airfoil2D.h"

namespace VM2D
{


	/*!
	\brief Класс, определяющий тип обтекаемого профиля

	Тип профиля:
	- профиль с криволинейными панелями.

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.11
	\date 07 августа 2022 г.
	*/

	class AirfoilCurv :
		public Airfoil
	{
	private: 
		///\brief Координаты центра панелей  
		std::vector<Point2D> rc_;

		/// Кривизна в центре панели
		std::vector<double> kc_;

		/// Производная от кривизны  в центре панели
		std::vector<double> dkc_;

		/// Вторая производная от кривизны  в центре панели
		std::vector<double> ddkc_;

		/// Кривизна в конце панели
		std::vector<double> k_;

		/// Производная от кривизны  в конце панели
		std::vector<double> dk_;

		/// Вторая производная от кривизны  в конце панели
		std::vector<double> ddk_;

	public:

		/// Конструктор
		AirfoilCurv(const World2D& W_, const size_t numberInPassport_)
			:Airfoil(W_, numberInPassport_)
		{ };

		AirfoilCurv(const Airfoil& afl) : Airfoil(afl) {};

		/// Деструктор
		virtual ~AirfoilCurv() { };

		/// Вычисление нормалей
		void CalcNrm();

		///Вычисляет габаритный прямоугольник профиля
		virtual void GetGabarits(double gap = 0.02) override;

		//далее -- реализация виртуальных функций
		virtual void ReadFromFile(const std::string& dir) override;
		virtual void Rotate(double alpha) override;
		virtual void Scale(const Point2D& factor) override;
		virtual void Move(const Point2D& dr) override;

		virtual bool IsPointInAirfoil(const Point2D& point) const override;

		virtual std::vector<double> getA(size_t p, size_t i, const Airfoil& otherAirfoil, size_t j) const override;
		virtual void calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const override;


		///\brief Возврат константной ссылки на координату центра панели
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getRc[size()] = getRc[0];		
		/// \n Это сделано для единообразия с векторами r_ и v_
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на координату центра панели
		const Point2D& getRc(size_t q) const
		{
			const size_t& sz = rc_.size();
			if (q < sz)
				return rc_[q];
			else
				return rc_[q-sz];
		};

		///\brief Возврат константной ссылки на кривизну в центре панели
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getKc[size()] = getKc[0];		
		/// \n Это сделано для единообразия с векторами r_ и v_
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на кривизну в центре панели 
		const double& getKc(size_t q) const
		{
			const size_t& sz = kc_.size();
			if (q < sz)
				return kc_[q];
			else
				return kc_[q - sz];
		};

		///\brief Возврат константной ссылки на производную кривизны в центре панели
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getDkc[size()] = getDkc[0];		
		/// \n Это сделано для единообразия с векторами r_ и v_
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на производную кривизны в центре панели 
		const double& getDkc(size_t q) const
		{
			const size_t& sz = dkc_.size();
			if (q < sz)
				return dkc_[q];
			else
				return dkc_[q - sz];			
		};

		///\brief Возврат константной ссылки на вторую производную кривизны в центре панели
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getDdkc[size()] = getDdkc[0];		
		/// \n Это сделано для единообразия с векторами r_ и v_
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на вторую производную кривизны в центре панели 
		const double& getDdkc(size_t q) const
		{		
			const size_t& sz = ddkc_.size();
			if (q < sz)
				return ddkc_[q];
			else
				return ddkc_[q - sz];
		};

		///\brief Возврат константной ссылки на кривизну в вершине профиля
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getK[size()] = getK[0];
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на кривизну в вершине профиля
		const double& getK(size_t q) const
		{
			const size_t& sz = k_.size();
			if (q < sz)
				return k_[q];
			else
				return k_[q - sz];			
		};

		///\brief Возврат константной ссылки на производную кривизны в вершине профиля
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getDk[size()] = getDk[0];
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на производную кривизны в вершине профиля
		const double& getDk(size_t q) const
		{			
			const size_t& sz = dk_.size();
			if (q < sz)
				return dk_[q];
			else
				return dk_[q - sz];
		};

		///\brief Возврат константной ссылки на вторую производную кривизны в вершине профиля
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getDdk[size()] = getDdk[0];		
		/// \n Это сделано для единообразия с векторами r_ и v_
		/// 
		/// \param[in] q номер панели профиля
		/// return константную ссылку на вторую производную кривизны в вершине профиля
		const double& getDdk(size_t q) const
		{
			const size_t& sz = ddk_.size();
			if (q < sz)
				return ddk_[q];
			else
				return ddk_[q - sz];
		};



		virtual void GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override;
#if defined(USE_CUDA)
		virtual void GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override;
#endif
		virtual void GetDiffVelocityI0I3ToWakeAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override {};

		virtual void GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const override;

		virtual void GetInfluenceFromSourcesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const override {};

		virtual void GetInfluenceFromSourceSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override {};

		virtual void GetInfluenceFromVortexSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override {};

        virtual void GetInfluenceFromVInfToPanel(std::vector<double>& vInfRhs) const override;
	};

} //namespace VM2D


#endif

/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DDeformable.h                                          |
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
\brief Заголовочный файл с описанием класса MechanicsDeformable
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef MECHANICSDEFORMABLE_H
#define MECHANICSDEFORMABLE_H

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;
	class Beam
	{
	public:
		double x0; //абсцисса начала балки
		double L;  //длина балки

		double rho, F, EJ;

		int R;
		
		//соб.частоты для единичной консольной балки
		const std::vector<double> unitLambda = { 1.8751, 4.6941, 7.8548, 10.9955, 14.1772, 17.2788, 20.4204 };
		std::vector<double> qCoeff;		

		//Интегралы от квадратов собственных форм
		const double intSqUnitShape = 0.25;

		//Собств.формы
		double shape(int n, double L, double x) const
		{
			double lam = unitLambda[n] / L;
			double C4 = (sin(unitLambda[n]) - sinh(unitLambda[n])) / (cos(unitLambda[n]) + cosh(unitLambda[n]));

			return 0.5 * (-cos(lam * (x - x0)) + cosh(lam * (x - x0)) - C4 * sin(lam * (x - x0)) + C4 * sinh(lam * (x - x0)));
		}

		std::vector<double> currentPhi, currentDPhi;

		const size_t nLastSteps = 10;
		std::vector<std::vector<double>> presLastSteps;

	public:
		Beam(double x0_, double L_, int R_) :
			rho(1000.0),
			F(0.02),
			EJ(1.11111),
			R(R_),
			x0(x0_),
			L(L_)
		{
			qCoeff.resize(R);
			currentPhi.resize(R);
			currentDPhi.resize(R);
			presLastSteps.reserve(nLastSteps);
		};

		double phi(int n, double t) const
		{
			return currentPhi[n];
		}

		void solveDU(int n, double dt)
		{
			double cm = rho * F;
			double ck = EJ * sqr(sqr(unitLambda[n] / L));
			double cq = qCoeff[n];

			double omega = sqrt(ck / cm);
			double cDamp = 0.005;

			double phiAst = currentPhi[n] + 0.5 * dt * currentDPhi[n];
			double psiAst = currentDPhi[n] + 0.5 * dt * (-ck * currentPhi[n] - cDamp * omega * currentDPhi[n] - cq) / cm;

			currentPhi[n] += dt * psiAst;
			currentDPhi[n] += dt * (-ck * phiAst - cDamp * omega * psiAst - cq) / cm;
		}


		double getTotalDisp(double x, double t) const //имитатор деформации упругой линии
		{ 
			double result = 0.0;
			for (int i = 0; i < R; ++i)
				result += phi(i, t) * shape(i, L, x);
			return result;
		}
	};


	struct ChordPanel //Хорда балки (упругая линия)
	{
		Point2D beg, end;					 //начало и конец
		std::pair<size_t, size_t> infPanels; //индексы панелей на обтекаемой поверхности балки под хордой и над хордой
		double rightSemiWidth;				 //полутолщина балки на правом конце
	};



	/*!
	\brief Класс, определяющий вид механической системы

	Деформируемое

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна

	\Version 1.12
	\date 14 января 2024 г.
	*/

	class MechanicsDeformable :
		public Mechanics
	{
	private:

		

	public:
		/// текущая скорость профиля
		Point2D& getVcm() { return Vcm; };		

		/// текущее отклонение профиля
		Point2D& getRcm() { return Rcm; };

		/// текущая угловая скорость профиля
		double& getWcm() { return Wcm; };

		/// текущий угол поворота профиля
		double& getPhicm() { return Phi; };

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
		MechanicsDeformable(const World2D& W_, size_t numberInPassport_);

		/// Деструктор
		~MechanicsDeformable() {};

		//далее -- реализации виртуальных функций
		virtual void GetHydroDynamForce() override;
		virtual Point2D VeloOfAirfoilRcm(double currTime) override;
		virtual Point2D PositionOfAirfoilRcm(double currTime) override;
		virtual double AngularVelocityOfAirfoil(double currTime) override;
		virtual double AngleOfAirfoil(double currTime) override;
		virtual void VeloOfAirfoilPanels(double currTime) override;
		virtual void ReadSpecificParametersFromDictionary() override;
		virtual void Move() override;

		/// хорда
		size_t indexOfUpperRightAngle;
		size_t indexOfUpperLeftAngle;
		size_t indexOfLowerRightAngle;
		size_t indexOfLowerLeftAngle;

		std::vector<ChordPanel> chord;    //геометрия упругой линии
		std::vector<double> upperShifts;  //превышения (по вертикали) точек верхней поверхности балки над хордой
		std::vector<double> lowerShifts;  //превышения (по вертикали) точек нижней поверхности балки над хордой

		std::unique_ptr<Beam> beam;       //сама балка

	};

}//namespace VM2D

#endif
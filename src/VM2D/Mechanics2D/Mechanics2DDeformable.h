/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\author Серебровская Екатерина Александровна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef MECHANICS2DDEFORMABLE_H
#define MECHANICS2DDEFORMABLE_H

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Вспомогательный класс Beam для описания упругой хорды деформируемой балки
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\author Серебровская Екатерина Александровна
	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class Beam
	{
	public:
		const bool fsi;
		
		const World2D& W;

		double x0; //абсцисса начала балки
		double L;  //длина балки

		double rho, F, EJ;

		int R;
		
		//соб.частоты для единичной консольной балки
		const std::vector<double> unitLambda = { 1.87510407, 4.69409113, 7.854757, 10.99554073, 14.1371683 };
		std::vector<double> qCoeff;

		//Интегралы от квадратов собственных форм
		const double intSqUnitShape = 0.25;

		//Собств.формы
		double shape(int n, double x) const
		{
			double lam = unitLambda[n] / L;
			double C4 = (sin(unitLambda[n]) - sinh(unitLambda[n])) / (cos(unitLambda[n]) + cosh(unitLambda[n]));

			return 0.5 * (-cos(lam * (x - x0)) + cosh(lam * (x - x0)) - C4 * sin(lam * (x - x0)) + C4 * sinh(lam * (x - x0)));
		}

		std::vector<double> currentPhi, currentDPhi;

		const size_t nLastSteps = 1;
		std::vector<std::vector<double>> presLastSteps;

	public:
		Beam(const World2D& W_, bool fsi_, double x0_, double L_, int R_);
		double phi(int n, double t) const;
		void solveDU(int n, double dt);
		void solveDU_RK(int n, double dt);
		double getTotalDisp(double x, double t) const;
		double getGivenLaw(double x, double t, double deformParam) const; //имитатор деформации упругой линии
	};


	struct ChordPanel //Хорда балки (упругая линия)
	{
		Point2D beg, end;					 //начало и конец
		std::pair<size_t, size_t> infPanels; //индексы панелей на обтекаемой поверхности балки под хордой и над хордой
		double rightSemiWidth;				 //полутолщина балки на правом конце
	};



	/*!
	\brief Класс, определяющий вид механической системы

	Деформируемое твердое тело

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\author Серебровская Екатерина Александровна
	\Version 1.14
	\date 6 марта 2026 г.
	*/

	class MechanicsDeformable :
		public Mechanics
	{
	private:
		double deformParam;
		bool fsi; //false for fish, true for turek
		

	public:
		/// текущая скорость профиля
		Point2D& getVcm() { return Vcm; };		

		/// текущее отклонение профиля
		Point2D& getRcm() { return Rcm; };

		/// текущая угловая скорость профиля
		double& getWcm() { return Wcm; };

		/// текущий угол поворота профиля
		double& getPhicm() { return Phi; };

		//bool getFsi() const { return fsi; }; //false for fish, true for turek

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

		std::vector<ChordPanel> initialChord;    //геометрия упругой линии в начальном состоянии
		std::vector<ChordPanel> chord;    //геометрия упругой линии
		std::vector<double> upperShifts;  //превышения (по вертикали) точек верхней половины торца балки над хордой
		std::vector<double> lowerShifts;  //превышения (по вертикали) точек нижней половины торца балки над хордой

		std::unique_ptr<Beam> beam;       //сама балка

		std::vector<std::vector<Point2D>> initialPossibleWays;

	};

}//namespace VM2D

#endif
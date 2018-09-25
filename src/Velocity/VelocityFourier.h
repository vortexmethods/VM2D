/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: VelocityFourier.h                                                |
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
\brief В РАБОТЕ Заголовочный файл с описанием класса VelocityFourier
\warning Пока в состоянии нерабочего прототипа
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#ifndef VELOCITYFOURIER_H
#define VELOCITYFOURIER_H

#include "unsupported/Eigen/FFT"

#include "defs.h"
#include "Velocity.h"

class World2D;

/*!
\brief В РАБОТЕ Класс, определяющий способ вычисления скоростей

\warning Пока в состоянии нерабочего прототипа

Способ вычисления скоростей
- при помощи решения вспомогательного уравнения для функции тока на грубой сетке, используя быстрое преобразование Фурье

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.2
\date 14 июня 2018 г.
*/
class VelocityFourier : public Velocity
{
private:
	numvector<int, 2> nNode;
	Point2D L;
	Point2D corner;

	numvector<int, 2> nCell;
	Point2D h;

	double(*Kernel)(double);

	std::vector<std::vector<size_t>> vic;  //списки вихрей, попавших в ячейки сетки
	
	size_t cellIj(int i, int j) const; //пересчет глобального номера ячейки по "локальным координатам" ячейки
	
	//double W(const Point2D& pt) const; //"размазанная" завихренность

	void FillVic();  //заполнение списков в массиве vic

	double w(int i, int j) const; //завихренность, приведенная к узлу (i,j)

	Point2D Skos(const Point2D& R, const Point2D& X) const
	{
		double dst2 = std::max(dist2(R, X), 1e-10);
		Point2D res = { -(R[1] - X[1]), (R[0] - X[0]) };
		res *= IDPI / dst2;
		return res;
	}

public:
	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу 	
	VelocityFourier(const World2D& W_) :
		Velocity(W_)
	{
		nNode = { 16, 16 };
		L = { 1.0, 1.0 };
		corner = { 0.0, 0.0 };

		Kernel = M4;

		nCell = { nNode[0] - 1, nNode[1] - 1 };
		h = { L[0] / nCell[0], L[1] / nCell[1] };
	};

	virtual ~VelocityFourier() {};

	/// \todo Реализовать 
	virtual void CalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius) {};
#if defined(USE_CUDA)
	virtual void GPUCalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius) {};
#endif

	virtual void CalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2) {};
#if defined(USE_CUDA)
	virtual void GPUCalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh = false) {};
#endif	
};




#endif
/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MeasureVP2D.h                                                    |
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
\brief Заголовочный файл с описанием класса MeasureVP
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.5   
\date 20 февраля 2019 г.
*/

#ifndef MEASUREVP_H
#define MEASUREVP_H

#include <memory>

#include "defs.h"

namespace VM2D
{
	class WakeDataBase;
	class World2D;

	/*!
	\brief Класс, отвечающий за вычисление поля скорости и давления в заданых точках для вывода

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.5
	\date 20 февраля 2019 г.
	*/

	class MeasureVP
	{
	private:
		/// Точки, которые считываются из файла (давление пишется в vtk-файлы)
		std::vector<Point2D> initialPoints;

		/// Точки, которые считываются из файла (давление пишется в vtk и csv-файлы)
		std::vector<Point2D> historyPoints;

		/// \brief Умный указатель на точки, в которых нужно вычислять в данный момент времени 
		/// Хранятся в виде "мнимой" системы вихрей, чтобы воспользоваться методом CalcConvVeloToSetOfPoints(...)
		std::unique_ptr<WakeDataBase> additionalWake;

		/// Скорости в нужных точках
		std::vector<Point2D> velocity;

		/// Радиусы вихрей в нужных точках
		std::vector<double> domainRadius;

		/// Давление в нужных точках
		std::vector<double> pressure;

		/// \todo Сделать учет давления на бесконечности
		//const double pInf = 1.0;

	protected:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

	public:
		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		MeasureVP(const World2D& W_);

		/// Деструктор
		~MeasureVP() { };

		/// \brief Чтение точек, в которых нужно посчитать давление и скорость
		///
		/// \param[in] dir константная ссылка на строку --- имя каталога, где лежит cчитываемый файл
		void ReadPointsFromFile(const std::string& dir);

		/// \brief Расчет и сохранение в файл поля скоростей и давления
		///
		/// \param[in] dir константная ссылка на строку, задающую каталог, куда сохранять файл с вихревым следом
		/// \param[in] step номер кадра для сохранения
		/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
		void CalcSaveVP(const std::string& dir, size_t step, timePeriod& time);

	};

}//namespace VM2D

#endif
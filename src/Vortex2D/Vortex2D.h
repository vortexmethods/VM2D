/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: Vortex2D.h                                                       |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.	                                                          |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса Vortex2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VORTEX2D_H_
#define VORTEX2D_H_

#include <iostream>

#include "Point2D.h"

/*!
\brief Класс, опеделяющий двумерный вихревой элемент
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Vortex2D
{
private:
	/// Радиус-вектор вихря
	Point2D pos;	
	
	/// Циркуляция вихря
	double gam;	

public:
	/// MPI-описатель типа
	static MPI_Datatype mpiVortex2D;

	/// Пустой конструктор
	Vortex2D() {};

	/// \brief Конструктор инициализации
	///
	/// \param[in] _r константная ссылка на радиус-вектор положения вихря
	/// \param[in] _g циркуляция (интенсивность) вихря
	Vortex2D(const Point2D& _r, const double _g)
		: pos(_r), gam(_g)	{ }
	
	/// Деструктор
	~Vortex2D() {};

	/// \brief Функция для доступа к радиус-вектору вихря
	/// \return ссылка на радиус-вектор вихря
	Point2D& r() { return pos; }

	/// \brief Функция для доступа для чтения к радиус-вектору вихря
	/// \return константная ссылка на радиус-вектор вихря
	const Point2D& r() const { return pos; }

	/// \brief Функция для доступа к циркуляции вихря
	/// \return ссылка на циркуляцию вихря
	double& g() { return gam; }

	/// \brief Функция для доступа для чтения к циркуляции вихря
	/// \return константная ссылка на циркуляцию вихря
	const double& g() const { return gam; }

	/// Cоздание MPI-описателя типа
	static void CreateMpiType();	
};

#endif
 
/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: VirtualWake.h                                                    |
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
\brief Заголовочный файл с описанием класса VirtualWake
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef VIRTUALWAKE_H
#define VIRTUALWAKE_H

//#include <string>
#include <vector>

//#include "defs.h"
#include "WakeDataBase.h"
//#include "Airfoil.h"
//#include "Vortex2D.h"

class World2D;
class Boundary;

/*!
\brief Класс, опеделяющий вихревой след (пелену)
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
class VirtualWake : public WakeDataBase
{
private:
	
public:
	/// Константная ссылка на решаемую задачу
	const World2D& W;

	/// Константная ссылка на границу, которой принадлежит данный виртуальный след
	const Boundary& bnd;
	
	/// Скорость вихрей виртуального следа конкретного профиля (равна Gamma/2) используется для расчета давления
	std::vector<Point2D> vecHalfGamma;

	/// Пара чисел: номер профиля и номер панели, на которой рожден виртуальный вихрь
	std::pair<size_t, size_t> aflPan;

	/// \brief Конструктор инициализации
	///
	/// \param[in] W_ константная ссылка на решаемую задачу	
	/// \param[in] bnd_ константная ссылка на границу	
	VirtualWake(const World2D& W_, const Boundary& bnd_)
		: W(W_), bnd(bnd_) { };
	
	/// Деструктор
	~VirtualWake(){ };

	/// \brief MPI-синхронизация вихревого следа
	///
	/// \todo В целях оптимизации можно подумать над .reserve()
	///
	/// Рассылка следа на все процессоры локальной группы процессоров, занятых решением данной задачи
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void WakeSynchronize();
		
	
	};

#endif

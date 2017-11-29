/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: Sheet.h                                                          |
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
\brief Заголовочный файл с описанием класса Sheet
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef SHEET_H
#define SHEET_H

#include <stddef.h>
#include <vector>


/*!
\brief Класс, опеделяющий слои на поверхности обтекаемого профиля

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Sheet
{
public:
	
	/// Список из характеристик свободного вихревого слоя на панелях
	std::vector<std::vector<double>> freeVortexSheet;

	/// Список из характеристик присоединенного вихревого слоя на панелях
	std::vector<std::vector<double>> attachedVortexSheet;

	/// Список из характеристик присоединенного слоя источников на панелях
	std::vector<std::vector<double>> attachedSourceSheet;

	/// Пустой конструктор 
	Sheet() { };

	/// Деструктор
	~Sheet() { };

	/// \brief Установка pазмерностей всех векторов и их обнуление
	///
	/// \param[in] np число панелей на профиле (внешняя размерность списков)
	/// \param[in] layerDim количество чисел, которыми характеризуются слои на каждой из панелей
	void SetLayersDim(size_t np, size_t layerDim);
};


#endif
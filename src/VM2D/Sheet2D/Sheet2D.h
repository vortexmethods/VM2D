/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Sheet2D.h                                                        |
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
\brief Заголовочный файл с описанием класса SheetV
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef SHEET_H
#define SHEET_H

#include <vector>
#include <cstring>
#include <memory>

namespace VM2D
{
	class World2D;

	/*!
	\brief Класс, опеделяющий слои на поверхности обтекаемого профиля

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	class Sheet
	{
	private:
		/// Список из характеристик свободного вихревого слоя на панелях
		std::vector<double> freeVortexSheet_;

		/// Список из характеристик присоединенного вихревого слоя на панелях
		std::vector<double> attachedVortexSheet_;

		/// Список из характеристик присоединенного слоя источников на панелях
		std::vector<double> attachedSourceSheet_;

	public:
		/// Константная ссылка на решаемую задачу
		const World2D& W;
		
		/// Число параметров в описании каждого слоя
		const int dim;

		/// \brief Конструктор 
		///
		/// \param[in] W_ константная ссылка на решаемую задачу	
		/// \param[in] dim_ число параметров в описании каждого слоя (1 - кусочно-постоянное решение, 2 - кусочно-линейное и т.п.)
		Sheet(const World2D& W_, int dim_);

		/// Деструктор
		~Sheet() { };

		/// \brief Установка pазмерностей всех векторов и их обнуление
		///
		/// \param[in] np число панелей на профиле (внешняя размерность списков)
		void SetLayers(size_t np);

		size_t getSheetSize() const
		{
			return freeVortexSheet_.size();
		}
		
		const double& freeVortexSheet(size_t n, size_t moment) const
		{
			return freeVortexSheet_[n*dim + moment];
		}
			
		const double& attachedVortexSheet(size_t n, size_t moment) const
		{
			return attachedVortexSheet_[n*dim + moment];
		}

		const double& attachedSourceSheet(size_t n, size_t moment) const
		{
			return attachedSourceSheet_[n*dim + moment];
		}

		double& freeVortexSheet(size_t n, size_t moment)
		{
			return freeVortexSheet_[n*dim + moment];
		}

		double& attachedVortexSheet(size_t n, size_t moment)
		{
			return attachedVortexSheet_[n*dim + moment];
		}

		double& attachedSourceSheet(size_t n, size_t moment)
		{
			return attachedSourceSheet_[n*dim + moment];
		}
	};

}//namespace VM2D

#endif
/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: paramsBH.h                                                       |
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
\brief Параметры метода Барнса - Хата для CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef PARAMSBH_H
#define PARAMSBH_H


namespace BH
{	
	
//#define OLD_OMP

	//Признак подсчета числа операций
	//Включаются ключами cmake
	//#define calcOp 

#define needTreeVrt

	/*!
	\brief Класс, содержащий параметры метода Баонса - Хата для CPU
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\Version 1.14
	\date 6 марта 2026 г.
	*/

	class params
	{
	public:
		/// Радиус вихревого элемента
		double eps;

		/// Квадрат радиуса вихревого элемента
		double eps2;

		/// Точность расчета скоростей:
		//  order = 1: монополь (0) + диполь (0)
		//  order = 2: монополь (0+1) + диполь (0) + квадруполь (0)
		//  order = 3: монополь (0+1+2) + диполь (0+1) + квадруполь (0) + октуполь (0)
		//  order = 4: монополь (0+1+2+3) + диполь (0+1+2) + квадруполь (0+1) + октуполь (0) + гексадекаполь (0)
		int order;

		/// Максимальное количество уровней дерева
		int NumOfLevelsVortex;

		/// Параметр точности 
		double theta;
	};
}//namespace BH

#endif
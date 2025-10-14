/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Params.h                                                         |
| Info: Source code of BH                                                     |
|                                                                             |
| This file is part of BH.                                                    |
| BH is free software: you can redistribute it and/or modify it               |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BH.  If not, see <http://www.gnu.org/licenses/>.                 |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовок класса, содержащего параметры решаемой задачи
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#ifndef PARAMS_H_
#define PARAMS_H_


namespace BH
{	
	
//#define OLD_OMP

	//Признак подсчета числа операций
	//Включаются ключами cmake
	//#define calcOp 

#define needTreeVrt

	/*!
	\brief Класс, содержащий параметры решаемой задачи
	\author Марчевский Илья Константинович
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\version 1.5
	\date 19 июня 2024 г.
	*/

	class params
	{
	public:
		/// Имя файла со списком вихрей
		///std::string vortexFile;

		/// Смещение вихревой пелены
		///Point2D wakeShift;

		/// Имя файла со списком вершин профиля
		///std::vector<std::string> airfoilFile;

		/// Имя файла со списком точек вычисления скорости
		///std::string vpFile;

		/// Радиус вихревого элемента
		double eps;
		/// Квадрат радиуса вихревого элемента
		double eps2;

		/// Скорость набегающего потока
		Point2D velInf;

		/// Название задачи для файла с результатом
		///std::string task;

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

		/// Сравнение с прямым методом
		///bool compare;

		/// Число повторных запусков кода
		///int runs;

		/// Чтение скоростей по прямому методу из файла (если его нет, он будет посчитан и сохранен)
		///bool BSfromFile;

		/// Сохранение скоростей по быстрому методу в файл
		///bool save;

		//Технические параметры
		/// Для номер уровня дерева, до которого генерируются OMP нити
		int maxLevelOmp;

		/// Признак останова для GMRES
		///double epsGMRES;

		/// Признак печати невязки на итерациях при решении интегрального уравнения
		///bool residualShow;


		/// 
		int outVortex = 1000;
		mutable int outVortexCoded; //вычисляется
		mutable int outCell; //вычисляется
		mutable std::vector<int> farCells; //вычисляется
	};
}//namespace BH

#endif
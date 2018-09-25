/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Parallel.h                                                       |
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
\brief Заголовочный файл с описанием класса Parallel и структуры parProp
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#ifndef PARALLEL_H
#define PARALLEL_H

#include <vector>

#include "mpi.h"

/*!
\brief Стрктура, содержащая параметры исполнения задачи в параллельном MPI-режиме

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/
struct parProp
{
	/// Список из чисел витков циклов, предназначенных для каждого процессора
	std::vector<int> len;
	
	/// Список, определяющий номер витка цикла, с которого должен начинать работу данный процессор
	std::vector<int> disp;

	/// Число витков, предназначенное текущему процессору
	int myLen;

	/// Индекс первого витка из числа витков, предназначенных текущему процессору
	int myDisp;

	/// Общее число витков, разделенное между всеми процессорами
	int totalLen;
};


/*!
\brief Класс, опеделяющий параметры исполнения задачи в параллельном MPI-режиме

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/
class Parallel
{
public:
	/// Коммуникатор для решения конкретной задачи
	MPI_Comm commWork;

	/// Локальный номер процессора, решающего конкретную задачу
	int myidWork;

	/// Число процессоров, решающих конкретную задачу
	int nProcWork;

	/// \brief Распределение задач по процессорам
	///
	/// \param[in] n число распределяемых витков цикла
	/// \param[in] bcastAll признак рассылки всей информации всем процессорам (по умолчанию false)
	/// \return структуру типа parProp, заполненную для текущего процессора
	parProp SplitMPI(size_t n, bool bcastAll = false) const;	
};
 
#endif

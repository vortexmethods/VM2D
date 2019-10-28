/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: WakeDataBase2D.h                                                 |
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
\brief Заголовочный файл с описанием класса WakeDataBase
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/

#ifndef WAKEDATABASE_H
#define WAKEDATABASE_H

#include "defs.h"

namespace VM2D
{
	class World2D;

	/*!
	\brief Класс, опеделяющий набор вихрей
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.6
	\date 28 октября 2019 г.
	*/

	class WakeDataBase
	{
	public:

		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// \brief Конструктор инициализации
		///
		/// \param[in] W_ константная ссылка на решаемую задачу	
		WakeDataBase(const World2D& W_)
			: W(W_) 
		{
#if defined(__CUDACC__) || defined(USE_CUDA)
			devNWake = 0;
			/// \todo Что имелось в виду? Зачем 2 раза?
			//devNWake = 0;
#endif
		};

		/// Список вихревых элементов
		std::vector<Vortex2D> vtx;

#if defined(USE_CUDA)

		/// Число вихрей, под хранение которых зарезервирована память на видеокарте
		mutable size_t devNWake;

		/// Адрес на видеокарте, где будет храниться коппия массива вихрей
		mutable double* devVtxPtr;

		/// Векторы для возврата вычисленных значений на хост и адреса для хранения результатов на видеокарте
		mutable std::vector<Point2D> tmpVel;
		mutable double* devVelPtr;

		mutable std::vector<double> tmpRad;
		mutable double* devRadPtr;

		mutable std::vector<double> tmpI0;
		mutable double* devI0Ptr;

		mutable std::vector<double> tmpI1;
		mutable double* devI1Ptr;

		mutable std::vector<Point2D> tmpI2;
		mutable double* devI2Ptr;

		mutable std::vector<Point2D> tmpI3;
		mutable double* devI3Ptr;

		//Данные для коллапса
		mutable std::vector<int> tmpNei;
		mutable int* devMeshPtr;
		mutable int* devNeiPtr;

		mutable std::vector<numvector<int, 2>> mesh;
#endif

		virtual ~WakeDataBase() {};

		/// \brief Считывание вихревого следа из файла 
		/// 
		/// \param[in] dir константная ссылка на строку, задающую каталог, где лежит файл с вихревым следом
		/// \param[in] fileName константная ссылка на строку, задающую имя файла с вихревым следом
		void ReadFromFile(const std::string& dir, const std::string& fileName);


		/// \brief MPI-синхронизация вихревого следа
		///
		/// Рассылка следа на все процессоры локальной группы процессоров, занятых решением данной задачи
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		void WakeSynchronize();
	};

}//namespace VM2D

#endif


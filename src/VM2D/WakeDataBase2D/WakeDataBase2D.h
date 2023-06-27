/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\brief Заголовочный файл с описанием класса WakeDataBase
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef WAKEDATABASE_H
#define WAKEDATABASE_H

#include "defs.h"
#include "Gpu2D.h"

#if defined(USE_CUDA)
	class VectorsForKnn;
#endif

namespace VM2D
{
	class World2D;

	/*!
	\brief Класс, опеделяющий набор вихрей
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
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
		{};

		/// Список вихревых элементов
		std::vector<Vortex2D/*, VM2D::MyAlloc<Vortex2D>*/> vtx;

#if defined(USE_CUDA)

		/// Адрес на видеокарте, где будет храниться копия массива вихрей
		mutable double* devVtxPtr;

		/// Векторы для возврата вычисленных значений на хост и адреса для хранения результатов на видеокарте
		mutable double* devVelPtr;
		mutable double* devRadPtr;

		mutable double* devI0Ptr;
		mutable float* devI0fPtr;
		mutable double* devI1Ptr;
		mutable double* devI2Ptr;
		mutable double* devI3Ptr;
		mutable float* devI3fPtr;

		//Данные для коллапса
		mutable int* devMeshPtr;
		mutable int* devNeiPtr;

		mutable std::vector<numvector<int, 2>> mesh;

		mutable VectorsForKnn* vecForKnn = nullptr;

#endif

		virtual ~WakeDataBase() {};

		/// \brief Считывание вихревого следа из файла 
		/// 
		/// \param[in] dir константная ссылка на строку, задающую каталог, где лежит файл с вихревым следом
		/// \param[in] fileName константная ссылка на строку, задающую имя файла с вихревым следом
		void ReadFromFile(const std::string& dir, const std::string& fileName);
			   
		/// \brief Сохранение вихревого следа в файл .vtk
		void SaveKadrVtk(const std::string& filePrefix = "Kadr") const;
	};

}//namespace VM2D

#endif


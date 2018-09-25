/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Wake.h                                                           |
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
\brief Заголовочный файл с описанием класса Wake
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#ifndef WAKE_H
#define WAKE_H

#include <string>
#include <vector>

#include "defs.h"
#include "Airfoil.h"
#include "Vortex2D.h"

class World2D;

/*!
\brief Класс, опеделяющий вихревой след (пелену)
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/
class Wake : public WakeDataBase
{
private:
	/// Вектор потенциальных соседей для будущего коллапса
	std::vector<int> neighb;

public:

	/// Константная ссылка на решаемую задачу
	const World2D& W;
	
	/// \brief Конструктор инициализации
	///
	/// \param[in] W_ константная ссылка на решаемую задачу	
	Wake(const World2D& W_)
		: W(W_) { };
	
	/// Деструктор
	~Wake(){ };

	/// \brief Считывание вихревого следа из файла 
	/// 
	/// \param[in] dir константная ссылка на строку, задающую каталог, где лежит файл с вихревым следом
	void ReadFromFile(const std::string& dir);
	
	/// \brief Сохранение вихревого следа в файл 
	/// 
	/// \param[in] outVtx константная ссылка на вектор вихрей, которые нужно вывести в файл
	/// \param[in] dir константная ссылка на строку, задающую каталог, куда сохранять файл с вихревым следом
	/// \param[in] step номер кадра для сохранения
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void SaveKadr(const std::vector<Vortex2D>& outVtx, const std::string& dir, size_t step, timePeriod& time) const;

	/// \brief Сохранение вихревого следа в файл .vtk
	/// 
	/// \param[in] dir константная ссылка на строку, задающую каталог, куда сохранять файл с вихревым следом
	/// \param[in] step номер кадра для сохранения
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void SaveKadrVtk(const std::vector<Vortex2D>& outVtx, const std::string& dir, size_t step, timePeriod& time) const;

	/// \brief MPI-синхронизация вихревого следа
	///
	/// \todo В целях оптимизации можно подумать над .reserve()
	///
	/// Рассылка следа на все процессоры локальной группы процессоров, занятых решением данной задачи
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void WakeSynchronize();

	/// \brief Проверка пересечения вихрями следа профиля при перемещении
	///
	/// Исполняется сразу для всех вихрей в пелене, осуществляет проверку для отдельного профиля
	/// \n Вихри, попавшие внутрь профиля, получают нулевую циркуляцию, а их "бывшая" циркуляция передается в вектор gammaThrough в структуру данных профиля
	/// 
	/// \param[in] newPos константная ссылка на вектор из новых положений вихрей в вихревом следе
	/// \param[in] isMoves признак того, что профиль подвижный
	/// \param[in] oldAfl константная ссылка контролируемый профиль до перемещения (используется, если у профиля стоит признак того, что он движется)
	/// \param[in,out] afl ссылка на контролируемый профиль (происходит изменение afl->gammaThrough)
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void Inside(const std::vector<Point2D>& newPos, Airfoil& afl, bool isMoves, const Airfoil& oldAfl);

	/// \brief Реструктуризация вихревого следа
	///
	/// Исполняется сразу для всех вихрей в пелене
	/// \n Вихри, находящиеся далеко от профилей, удаляются
	/// \n Вихри, которые сильно сблизились, коллапсируются
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void Restruct(timePeriod& time);

	/// \brief Зануление далеко улетевших вихрей
	/// \return число вихрей в дальнем следе, которые занулены
	int KillFar();

	/// \brief Исключение нулевых и мелких вихрей
	/// \return число исключенных вихрей
	size_t RemoveZero();

	/// \brief Проверка проникновения точки через границу профиля
	/// 
	/// \param[in] newPos константная ссылка на смещенное (новое) положение
	/// \param[in] oldPos константная ссылка на несмещенное (старое) положение
	/// \param[in] afl константная ссылка на контролируемый профиль
	/// \param[out] panThrough номер "протыкаемой" панели
	/// return признак пересечения профиля
	bool MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, size_t& panThrough);

	/// \brief Проверка проникновения точки через границу профиля
	/// 
	/// \param[in] newPos константная ссылка на смещенное (новое) положение вихря
	/// \param[in] oldPos константная ссылка на несмещенное (старое) положение вихря
	/// \param[in] oldAfl константная ссылка на состояние контролируемого профиля до перемещения
	/// \param[in] afl константная ссылка на контролируемый профиль
	/// \param[out] panThrough номер "протыкаемой" панели
	/// return признак пересечения профиля
	bool MoveInsideMovingBoundary(const Point2D& newPos, const Point2D& oldPos, const Airfoil& oldAfl, const Airfoil& afl, size_t& panThrough);

	/// \brief Поиск ближайшего соседа
	/// \param[in] type тип коллапса: 
	/// - 0 --- без приоритета знаков
	/// - 1 --- коллапсировать только вихри разных знаков
	/// - 2 --- коллапсировать только вихри одного знака
	void GetPairs(int type);

#if defined(USE_CUDA)
	void GPUGetPairs(int type);
#endif

	/// \brief Коллапс вихрей
	/// \param[in] type тип коллапса: 
	/// - 0 --- без приоритета знаков
	/// - 1 --- коллапсировать только вихри разных знаков
	/// - 2 --- коллапсировать только вихри одного знака
	/// \param[in] times число проходов алгоритма коллапса
	/// \return число зануленных вихрей
	int Collaps(int type, int times);
	
	};

#endif

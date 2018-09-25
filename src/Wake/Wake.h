/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
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
\version 1.1
\date 2 апреля 2018 г.
*/

#ifndef WAKE_H
#define WAKE_H

#include "Airfoil.h"
#include "gpudefs.h"


#include <array>

class Boundary;
class gpu;

/*!
\brief Класс, опеделяющий вихревой след (пелену)
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/
class Wake
{
private:
	/// Вектор потенциальных соседей для будущего коллапса
	std::vector<int> neighb;

public:
	/// Список вихревых элементов
	std::vector<Vortex2D> vtx;
	IFCUDA(mutable double* devWakePtr);
	IFCUDA(mutable size_t devNWake);

	/// Константная ссылка на параметры дискретизации вихревого следа
	const WakeDiscretizationProperties& param;
	
	/// Константная ссылка на параметры исполнения в параллельном режиме
	const Parallel& parallel;

	/// Ссылка на объект, управляющий графическим ускорителем
	gpu& cuda;

	/// Константная ссылка на список умных казателей на обтекаемые профили
	const std::vector<std::unique_ptr<Airfoil>>& airfoils;
	
	/// Константная ссылка на список умных казателей на граничные условия
	const std::vector<std::unique_ptr<Boundary>>& boundary;

	/// \brief Конструктор инициализации
	///
	/// \param[in] param_ константная ссылка на параметры дискретизации вихревого следа
	/// \param[in] parallel_ константная ссылка на параметры исполнения в параллельном режиме
	/// \param[in] cuda_ ссылка на объект, управляющий графическим ускорителем
	/// \param[in] airfoils_ константная ссылка на вектор указателей на профили
	/// \param[in] boundary_ константная ссылка на вектор указателей на грфничные условия
	Wake(const WakeDiscretizationProperties& param_, const Parallel& parallel_, gpu& cuda_, const std::vector<std::unique_ptr<Airfoil>>& airfoils_, const std::vector<std::unique_ptr<Boundary>>& boundary_)
		: param(param_), parallel(parallel_), cuda(cuda_), airfoils(airfoils_), boundary(boundary_) { };
	
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
	void SaveKadrVtk(const std::string& dir, size_t step, timePeriod& time) const;

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
	/// \param[in,out] afl ссылка на контролируемый профиль (происходит изменение afl->gammaThrough)
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void Inside(const std::vector<Point2D>& newPos, Airfoil& afl);

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

	/// \brief Коллапс вихрей
	/// \param[in] type тип коллапса: 
	/// - 0 --- без приоритета знаков
	/// - 1 --- коллапсировать только вихри разных знаков
	/// - 2 --- коллапсировать только вихри одного знака
	/// \param[in] times число проходов алгоритма коллапса
	/// \return число зануленных вихрей
	int Collaps(int type, int times);


	/// \brief Проверка пересечения вихрями следа профиля при подвижном профиле
	///
	/// Исполняется сразу для всех вихрей в пелене, осуществляет проверку для отдельного профиля
	/// \n Вихри, попавшие внутрь профиля, получают нулевую циркуляцию, а их "бывшая" циркуляция передается в вектор gammaThrough в структуру данных профиля
	/// 
	/// \param[in] newPos константная ссылка на вектор из новых положений вихрей в вихревом следе
	/// \param[in] oldAfl константная ссылка контролируемый профиль до перемещения
	/// \param[in,out] afl ссылка на контролируемый профиль (происходит изменение afl->gammaThrough)
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void InsideMovingBoundary(const std::vector<Point2D>& newPos, const Airfoil& oldAfl, Airfoil& afl);

	};

#endif

/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary.h                                                       |
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
\brief Заголовочный файл с описанием класса Boundary
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include "Eigen/Dense"

#include "Airfoil.h"
#include "Point2D.h"
#include "Sheet.h"
#include "Vortex2D.h"
#include "VirtualWake.h"

class World2D;

/*!
\brief Абстрактный класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

class Boundary
{
protected:
	/// Константная ссылка на решаемую задачу
	const World2D& W;

	/// Номер профиля в паспорте
	const size_t numberInPassport;	

public:	

	/// Константная ссылка на профиль 
	/// \n инициализируется автоматом в конструкторе
	const Airfoil& afl;

	/// Константная ссылка на вектор из начал (и концов) панелей
	const std::vector<Point2D>& CC;

	/// Число вихрей, рождаемых на каждой панели профиля и формирующих виртуальный вихревой след
	int vortexPerPanel;
	
	/// Виртуальный вихревой след конкретного профиля
	VirtualWake virtualWake;

	/// \brief MPI-синхронизация виртуального вихревого следа
	///
	/// \todo В целях оптимизации можно подумать над .reserve()
	///
	/// Рассылка следа на все процессоры локальной группы процессоров, занятых решением данной задачи
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	void VirtualWakeSynchronize();

	/// \brief Размерность параметров каждого из слоев на каждой из панелей
	///
	/// Указывает, сколькими числами задается интенсивность каждого из слоев на каждой панели:
	/// - 1 --- одно число --- задается только среднее значение;
	/// - 2 --- два числа --- задается среднее значение и "наклон".
	size_t sheetDim;

	/// Слои на профиле
	Sheet sheets;

	/// Слои на профиле с предыдущего шага
	Sheet oldSheets;

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи
	/// \param[in] sheetDim_ размерность параметра слоя на каждой панели профиля (сколько чисел используется, чтобы описать слой на каждой панели);
	Boundary(const World2D& W_, size_t numberInPassport_, int sheetDim_);
	
	/// Деструктор
	virtual ~Boundary() { };
	
	/// \brief Генерация блока матрицы
	///
	/// Генерирует следующие компоненты матрицы:
	/// - диагональный блок матрицы --- влияние данного профиля на самого себя;
	/// - нижнюю строку для матрицы для данного профиля;
	/// - правый столбец матрицы для данного профиля.
	/// 
	/// \param[out] matr ссылка на генерируемую матрицу
	/// \param[out] lastLine ссылка на нижнюю строку
	/// \param[out] lactCol ссылка на правый столбец
	virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol) = 0;



	/// \brief Генерация блока матрицы влияния от другого профиля того же типа
	///
	/// \todo ДОДЕЛАТЬ ОПИСАНИЕ
	/// \todo Пока считается, что граничные условия одинаковые
	virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr) = 0;


	/// \brief Генерация вектора влияния вихревого следа на профиль
	///
	/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
	/// 
	/// \param[out] wakeVelo ссылка на вектор влияния вихревого следа на профиль
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const = 0;
#if defined(USE_CUDA)
	virtual void GPUGetWakeInfluence(std::vector<double>& wakeVelo) const = 0;
#endif

	/// \brief Вычисление конвективных скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как от слоев
	///
	/// Вычисляет конвективные скорости в наборе точек, которые вызваны влиянием завихренности и источников на профиле как от слоев
	/// 
	/// \param[in] points константная ссылка на набор точек, в которых вычисляются скорости
	/// \param[out] velo ссылка на вектор скоростей, которые приобретают точки из-за влияния завихренности и источников на профиле
	/// 
	/// \warning velo --- накапливается!
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const = 0;

	/// \brief Вычисление конвективных скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как виртуальных вихрей
	///
	/// Вычисляет конвективные скорости в наборе точек, которые вызваны влиянием завихренности и источников на профиле как виртуальных вихрей
	/// 
	/// \param[in] pointsDb константная ссылка на базу данных вихрей, в точках которых вычисляются скорости
	/// \param[out] velo ссылка на вектор скоростей, которые приобретают точки из-за влияния завихренности и источников на профиле
	/// 
	/// \warning velo --- накапливается!
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const = 0;
#if defined(USE_CUDA)
	virtual void GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const = 0;
#endif

	/// \brief Заполнение в правой части влияния набегающего потока, следа и присоединенных слоев, действующих от самого себя
	///
	/// Заполняет блок правой части матрицы СЛАУ, обеспечивающей удовлетворение граничного условия, соответствующий данному профилю.
	/// Учитывается влияние набегающего потока, следа и присоединенных слоев, действующих от самого себя
	/// 	
	/// \param[in] V0 константная ссылка на вектор набегающего потока
	/// \param[out] rhs ссылка на блок вектора правой части
	/// \param[out] lastRhs указатель на последний элемент вектора правой части
	/// \param[in] move является ли профиль подвижным
	/// \param[in] deform является ли профиль деформируемым
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs, bool move, bool deform) = 0;

	/// \brief Заполнение в правой части влияния присоединенных слоев, действующих на один профиль от другого
	///
	/// Заполняет блок правой части матрицы СЛАУ, обеспечивающий удовлетворение граничного условия, .
	/// Учитывается влияние присоединенных слоев от другого профиля
	/// 	
	/// \param[out] otherAirfoil константная ссылка профиль, вклад от которого в правую часть требуется вычислить
	/// \param[out] rhs ссылка на блок вектора правой части, соответствующего профилю, на который учитывается влияние
	///
	/// \todo Пока считается, что граничные условия одинаковые
	/// \warning Меняет rhs путем прибавления
	virtual void FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs) = 0;


	/// \brief Возврат размерности вектора решения 
	///
	/// (без учета регуляризирующей переменной)
	///
	/// \return размерность вектора решения
	virtual size_t GetUnknownsSize() const = 0;

	/// \brief Пересчет решения на интенсивность вихревого слоя и на рождаемые вихри на конкретном профиле
	///
	/// 1) Приводит решение к интенсивности вихревого слоя и записывает его в sheets.freeVortexSheet:
	///
	/// - если неизвестное --- интенсивность вихря, то он "размазывается" по панели;
	/// - если неизвестное --- интенсивность слоя, то она передается непостредственно.
	///
	/// 2) Приводит интенсивность вихревого слоя к рождаемым вихрям, а также вычисляет их положения
	/// \param[in] sol вектор решения СЛАУ
	virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol) = 0;

	/// \brief Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
	virtual void ComputeAttachedSheetsIntensity() = 0;
};

#endif

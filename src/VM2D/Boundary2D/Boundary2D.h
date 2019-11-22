/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary2D.h                                                     |
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
\brief Заголовочный файл с описанием класса Boundary
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <memory>
#include "Sheet2D.h"
#include "VirtualWake2D.h"

namespace VM2D
{

	class World2D;
	class Airfoil;

	/*!
	\brief Абстрактный класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.7
	\date 22 ноября 2019 г.
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

		/// Минимальное число вихрей, рождаемых на панели профиля и формирующих виртуальный вихревой след 
		int minVortexPerPanel;

		/// Число вихрей, рождаемых на каждой панели профиля (формируется после решения СЛАУ)
		std::vector<std::pair<std::vector<Vortex2D>::iterator, std::vector<Vortex2D>::iterator>> vortexBeginEnd;

		/// Виртуальный вихревой след конкретного профиля
		VirtualWake virtualWake;

		/// \brief MPI-синхронизация виртуального вихревого следа
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
		/// \param[in] currentRow номер строки (или столбца), с которого начинается текущий блок в общей матрице
		virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol) = 0;

		/// \brief Генерация блока матрицы, состоящей из интегралов от (r-xi)/|r-xi|^2, влияние профиля самого на себя
		///
		/// \param[out] IQ ссылка на генерируемую матрицу
		/// \param[in] currentRow номер строки (или столбца), с которого начинается текущий блок в общей матрице
		virtual void FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) = 0;

		/// \brief Генерация блока матрицы влияния от другого профиля того же типа
		///
		/// Генерирует блок матрицы влияния от другого профиля того же типа
		/// 
		/// \param[in] otherBoundary константная ссылка на граничное условие на втором профиле
		/// \param[out] matr ссылка на генерируемый блок матрицы
		/// \param[in] currentRow номер строки, с которого начинается текущий блок в общей матрице
		/// \param[in] currentCol номер столбца, с которого начинается текущий блок в общей матрице
		/// \todo Пока считается, что граничные условия одинаковые
		virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr) = 0;

		/// \brief Генерация блока матрицы, состоящей из интегралов от (r-xi)/|r-xi|^2, влияние одного профиля на другой
		///
		/// Генерирует блок матрицы влияния от другого профиля того же типа
		/// 
		/// \param[in] otherBoundary константная ссылка на граничное условие на втором профиле
		/// \param[out] matr ссылка на генерируемый блок матрицы
		/// \param[in] currentRow номер строки, с которого начинается текущий блок в общей матрице
		/// \param[in] currentCol номер столбца, с которого начинается текущий блок в общей матрице
		/// \todo Пока считается, что граничные условия одинаковые
		virtual void FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) = 0;

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

		/// \brief Вычисление конвективной скорости в точках расположения виртуальных вихрей
		///
		/// Вычисление производится в точках расположения виртуальных вихрей.
		/// Скорости находятся как скорости соответствующих точек профиля плюс половина интенсивности 
		/// присоединенного вихревого слоя, умноженная на касательную к профилю на данной панели.
		///
		/// \param[out] velo ссылка на заполняемый список скоростей
		/// \warning Массив velo заполняется путем присвоения, а не прибавления значений
		virtual void GetConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const = 0;

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
		size_t GetUnknownsSize() const;


		/// \brief Пересчет решения на интенсивность вихревого слоя и на рождаемые вихри на конкретном профиле
		///
		/// 1) Приводит решение к интенсивности вихревого слоя и записывает его в sheets.freeVortexSheet:
		///
		/// - если неизвестное --- интенсивность вихря, то он "размазывается" по панели;
		/// - если неизвестное --- интенсивность слоя, то она передается непосредственно.
		///
		/// 2) Приводит интенсивность вихревого слоя к рождаемым вихрям, а также вычисляет их положения
		/// \param[in] sol вектор решения СЛАУ
		virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol) = 0;

		/// \brief Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
		virtual void ComputeAttachedSheetsIntensity() = 0;


		/// \brief Вычисление влияния части подряд идущих вихрей из вихревого следа на прямолинейную панель для правой части
		///
		/// Вычисляет влияния части подряд идущих вихрей из вихревого следа на прямолинейную панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона вихрей
		/// \param[in] count длина диапазона вихрей
		/// \param[out] wakeRhs ссылка на вектор полученных влияние для правой части СЛАУ
		virtual void GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const = 0;
	
		/// \brief Вычисление влияния части подряд источников из области течения на прямолинейную панель для правой части
		///
		/// Вычисляет влияния части подряд идущих источников из области течения на прямолинейную панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона источников
		/// \param[in] count длина диапазона источников
		/// \param[out] wakeRhs ссылка на вектор полученных влияние для правой части СЛАУ
		virtual void GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const = 0;


		/// \brief Вычисление влияния части подряд идущих вихрей из вихревого следа на криволинейную панель для правой части
		///
		/// Вычисляет влияния части подряд идущих вихрей из вихревого следа на криволинейную панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона вихрей
		/// \param[in] count длина диапазона вихрей
		/// \param[out] wakeRhs ссылка на вектор полученных влияние для правой части СЛАУ
		virtual void GetInfluenceFromVorticesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const = 0;

		/// \brief Вычисление влияния части подряд идущих источников из области течения на криволинейную панель для правой части
		///
		/// Вычисляет влияния части подряд идущих источников из области течения на криволинейную панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона источников
		/// \param[in] count длина диапазона источников
		/// \param[out] wakeRhs ссылка на вектор полученных влияние для правой части СЛАУ
		virtual void GetInfluenceFromSourcesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const = 0;


};

}//namespace VM2D

#endif

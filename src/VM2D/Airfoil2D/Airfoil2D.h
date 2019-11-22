/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Airfoil2D.h                                                      |
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
\brief Заголовочный файл с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include "defs.h"

namespace VM2D
{

	class World2D;
	class WakeDataBase;

	/*!
	\brief Абстрактный класс, определяющий обтекаемый профиль
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class Airfoil
	{
	protected:
		/// Координаты начал панелей  
		std::vector<Point2D> r_;

		/// Скорости начал панелей  
		std::vector<Point2D> v_;

	public:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// Номер профиля в паспорте
		const size_t numberInPassport;

		/// Положение центра масс профиля
		Point2D rcm;

		/// Поворот профиля
		double phiAfl;

		/// Масса профиля
		double m;

		/// Полярный момент инерции профиля относительно центра масс
		double J;

		///\brief Возврат константной ссылки на вершину профиля  
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getR[size()] = getR[0];		
		/// \n Это позволяет удобно обращаться к getR(i) и getR(i+1) как к началу и концу i-й панели
		/// 
		/// \param[in] q номер вершины профиля
		/// return константную ссылку на вершину профиля  
		const Point2D& getR(size_t q) const
		{
			return (q < r_.size()) ? r_[q] : r_[0];
		};

		///\brief Возврат константной ссылки на скорость вершины профиля  
		///
		/// Организовано "зацикливание" в сторону увеличения индекса, т.е. getV[size()] = getV[0];		
		/// \n Это позволяет удобно обращаться к getV(i) и getV(i+1) как к скоростям начала и конца i-й панели
		/// 
		/// \param[in] q номер вершины профиля
		/// return константную ссылку на скорость вершины профиля  
		const Point2D& getV(size_t q) const
		{
			return (q < v_.size()) ? v_[q] : v_[0];
		};

		/// \brief Установка постоянной скорости всех вершин профиля
		///
		/// \param[in] vel константная ссылка на величину устанавливаемой скорости
		void setV(const Point2D& vel)
		{
			v_.clear();
			v_.resize(r_.size(), vel);
		}

		/// \brief Установка скоростей всех вершин профиля
		///
		/// \param[in] vel константная ссылка на вектор величин скоростей вершин профиля
		void setV(const std::vector<Point2D>& vel)
		{
			v_.clear();
			v_.insert(v_.end(), vel.begin(), vel.end());				
		}

		/// \brief Возврат количества панелей на профиле
		///
		/// \return количество панелей на профиле
		size_t getNumberOfPanels() const { return r_.size(); };
		
		/// Признак разворота нормалей (для расчета внутренних течений)
		bool inverse;

		/// Указатель на девайсе, где хранятся вершины профиля
		IFCUDA(mutable double* devRPtr);

		/// Указатель на девайсе, где хранится правая часть матрицы
		IFCUDA(mutable double* devRhsPtr);

		/// Указатель на хосте, где хранится временная часть матрицы, полученная с девайса
		IFCUDA(mutable std::vector<double> tmpRhs);

		/// Указатель на девайсе, где хранятся интенсивности свободного вихревого слоя на панелях
		IFCUDA(mutable double* devFreeVortexSheetPtr);

		/// Указатель на девайсе, где хранятся интенсивности присоединенного вихревого слоя на панелях
		IFCUDA(mutable double* devAttachedVortexSheetPtr);

		/// Указатель на девайсе, где хранятся интенсивности присоединенного слоя источников на панелях
		IFCUDA(mutable double* devAttachedSourceSheetPtr);

		/// Указатель на девайсе, где хранится вектор (по панелям) для силы вязкого трения
		IFCUDA(mutable double* devViscousStressesPtr);
			   
		/// Указатель на хосте, где хранится временная часть вектора (по панелям) для силы вязкого трения
		IFCUDA(mutable std::vector<double> tmpViscousStresses);

		/// \brief Нормали к панелям профиля
		///
		/// Нормали задаются внешними, нормированными на единицу
		std::vector<Point2D> nrm;

		/// \brief Касательные к панелям профиля
		///
		/// Касательные соответствуют обходу профиля против часовой стрелки, задаются нормированными на единицу
		std::vector<Point2D> tau;

		/// Длины панелей профиля
		std::vector<double> len;

		/// Касательные напряжения на панелях профиля
		std::vector<double> viscousStress;

		Point2D lowLeft; ///< Левый нижний угол габаритного прямоугольника профиля
		Point2D upRight; ///< Правый верхний угол габаритного прямоугольника профиля

		///\brief Суммарные циркуляции вихрей, пересекших панели профиля на прошлом шаге 
		///
		/// Используются в правой части системы, чтобы компенсировать вихри, проникшие в профиль
		std::vector<double> gammaThrough;

		/// Конструктор
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи
		Airfoil(const World2D& W_, const size_t numberInPassport_);


		/// Деструктор
		virtual ~Airfoil() { };


		/// \brief Проверка, идет ли вершина i следом за вершиной j
		///
		/// \param[in] i проверяемая вершина
		/// \param[in] j контрольная вершина
		/// \return true, если i-я вершина следует зп j-й в порядке обхода профиля
		bool isAfter(size_t i, size_t j) const;

		/// \brief Поворот профиля 
		///
		/// Поворачивает профиль на угол \f$ \alpha \f$ вокруг центра масс
		/// \param[in] alpha угол поворота против часовой стрелки в радианах
		virtual void Rotate(double alpha) = 0;

		/// \brief Масштабирование профиля 
		///
		/// Масштабирует профиль в factor раз относительно центра масс
		/// \param[in] factor масштабный коэффициент
		virtual void Scale(double factor) = 0;

		/// \brief Перемещение профиля 
		///
		/// Плоскопараллельно перемещает профиль на вектор \f$ \overrightarrow{dr} \f$
		///
		/// \param[in] dr константная ссылка на вектор перемещения
		virtual void Move(const Point2D& dr) = 0;

		/// \brief Вычисляет габаритный прямоугольник профиля
		///
		/// Заполняет значения полей lowLeft и upRight по габаритному прямоугольнику профиля с учетом зазора
		///
		/// \param[in] gap относительная величина зазора в долях от размера габаритного прямоугольника (по умолчанию 0.02, что соответствует 2 %)
		virtual void GetGabarits(double gap = 0.02) = 0;

		/// \brief Определяет, находится ли точка с радиус-вектором \f$ \vec r \f$ внутри габаритного прямоугольника профиля
		///
		/// \param[in] r константная ссылка на текущее положение точки
		/// \return true, если точка внутри габаритного прямоугольника
		bool isInsideGabarits(const Point2D& r) const;

		/// \brief Определяет, находится ли точка с радиус-вектором \f$ \vec r \f$ вне габаритного прямоугольника профиля
		///
		/// \param[in] r константная ссылка на текущее положение точки
		/// \return true, если точка вне габаритного прямоугольника
		bool isOutsideGabarits(const Point2D& r) const;

		/// \brief Определяет, находится ли точка с радиус-вектором \f$ \vec r \f$ внутри профиля
		///
		/// \param[in] point константная ссылка на текущее положение точки
		/// \return true, если точка внутри профиля
		virtual bool IsPointInAirfoil(const Point2D& point) const = 0;

		/// \brief Считывание профиля из файла
		///
		/// Считывание геометрии профиля из файла, вычисление всех прочих параметров профиля
		/// \n После загрузки из файла профиль поворачивается на нужный угол и масштабируется на нужный коэффициент
		/// \warning Сейчас масса, момент инерции и скорости вершин зануляются.
		///
		/// \param[in] dir константная ссылка на строку --- имя каталога, где лежит cчитываемый файл
		virtual void ReadFromFile(const std::string& dir) = 0;

		/// \brief Вычисление коэффициентов матрицы A для расчета влияния панели на панель
		///
		/// \param[in] p размерность матрицы - результата
		/// \param[in] i номер панели, на которую оказывается влияние
		/// \param[in] otherAirfoil константная ссылка на профиль, от которого рассчитывается влияние
		/// \param[in] j номер влияющей панели
		/// return соответствующий блок матрицы СЛАУ, вытянутый в линию 
		virtual std::vector<double> getA(size_t p, size_t i, const Airfoil& otherAirfoil, size_t j) const = 0;

		/// \brief Вычисление коэффициентов матрицы, состоящей из интегралов от (r-xi)/|r-xi|^2 
		///
		/// \param[in] p размерность матрицы - результата
		/// \param[in] otherAirfoil константная ссылка на профиль, от которого рассчитывается влияние
		/// return соответствующий блок матрицы СЛАУ, вытянутый в линию 
		virtual void calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const = 0;

		/// \brief Вычисление влияния присоединенных слоев от другого профиля (константные базисные функции)
		///
		/// Расчет интегралов в правой части с константными базисными функциями
		///
		/// \param[out] wakeVelo массив значений влияний для панелей профиля
		virtual void getInfAttFromOther0(std::vector<double>& attOtherVelo, const Airfoil& otherAirfoil, size_t currentRow, size_t currentCol) const = 0;

		/// \brief Вычисление влияния присоединенных слоев от другого профиля (линейные базисные функции)
		///
		/// Расчет интегралов с константными и линейными базисными функциями
		///
		/// \param[out] wakeVelo массив значений влияний для панелей профиля
		virtual void getInfAttFromOther1(std::vector<double>& attOtherVelo, const Airfoil& otherAirfoil, size_t currentRow, size_t currentCol) const = 0;



		/// \brief Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
		///
		/// Вычисляет диффузионные скорости в наборе точек, которые обусловленных геометрией профиля, и вычисляет вязкое трение
		/// 
		/// \param[in] pointsDb константная ссылка на базу данных вихрей, в которых вычисляются скорости
		/// \param[in] domainRadius ссылка на радиусы вихрей
		/// \param[out] I0 ссылка на вектор знаменателей диффузионных скоростей, которые приобретают точки из-за влияния геометрии профиля
		/// \param[out] I3 ссылка на вектор числителей диффузионных скоростей, которые приобретают точки из-за влияния геометрии профиля
		/// 
		/// \warning Векторы I0, I3 --- накапливаются!
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		virtual void GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) = 0;
#if defined(USE_CUDA)
		virtual void GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) = 0;
#endif

		/// \brief Вычисление влияния части подряд идущих вихрей из вихревого следа на панель для правой части
		///
		/// Вычисляет влияния части подряд идущих вихрей из вихревого следа на панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона вихрей
		/// \param[in] count длина диапазона вихрей
		/// \param[out] panelRhs ссылка на вектор полученного влияния для правой части СЛАУ для конкретной панели
		virtual void GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const = 0;
				
	
		/// \brief Вычисление влияния части подряд идущих источников из области течения на панель для правой части
		///
		/// Вычисляет влияния части подряд идущих источников из области течения на панель для правой части
		/// 
		/// \param[in] panel номер панели профиля, на которую считается влияние
		/// \param[in] ptr указатель на начало диапазона источников
		/// \param[in] count длина диапазона источников
		/// \param[out] panelRhs ссылка на вектор полученного влияния для правой части СЛАУ для конкретной панели
		virtual void GetInfluenceFromSourcesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const = 0;

};

}//namespace VM2D

#endif

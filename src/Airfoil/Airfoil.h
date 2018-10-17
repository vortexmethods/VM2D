/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Airfoil.h                                                        |
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
\brief Заголовочный файл с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <string>
#include <vector>

#include "Point2D.h"
#include "WakeDataBase.h"

class World2D;

/*!
\brief Абстрактный класс, определяющий обтекаемый профиль
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
class Airfoil
{
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
	
	/// Число панелей на профиле
	size_t np;    
	
    ///\brief Координаты начал панелей  
	///
	/// Всегда обеспечивается условие r[np] = r[0];
	/// \n вектор r содержит (np+1) элемент, где np --- число панелей.
	/// \n Это позволяет удобно обращаться к r[i] и r[i+1] как к началу и концу i-й панели
	std::vector<Point2D> r;
	IFCUDA(mutable double* devRPtr);
	IFCUDA(mutable double* devRhsPtr);
	IFCUDA(mutable std::vector<double> tmpRhs);

	IFCUDA(mutable double* devFreeVortexSheetPtr);
	IFCUDA(mutable double* devAttachedVortexSheetPtr);
	IFCUDA(mutable double* devAttachedSourceSheetPtr);
	
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
	
	///\brief Скорости начал панелей  
	///
	/// Всегда обеспечивается условие v[np] = v[0];
	/// \n вектор v содержит (np+1) элемент, где np --- число панелей.
	/// \n Это позволяет удобно обращаться к v[i] и v[i+1] как к скорости начала и конца i-й панели
	/// \warning при чтении из файла зашивается нулями
	std::vector<Point2D> v;

	///\brief Суммарные циркуляции вихрей, пересекших панели профиля на прошлом шаге 
	///
	/// Используются в правой части системы, чтобы компенсировать вихри, проникшие в профиль
	std::vector<double> gammaThrough;
	
	/// Конструктор
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи
	Airfoil(const World2D& W_, const size_t numberInPassport_);

	/// Конструктор копирования
	/// \todo откомментировать	
	Airfoil(const Airfoil& afl) : np(0), W(afl.W), numberInPassport(afl.numberInPassport)
	{	
		rcm = afl.rcm;
		np = afl.np;
		phiAfl = afl.phiAfl;
		r.insert(r.begin(), afl.r.begin(), afl.r.end());
		nrm.insert(nrm.begin(), afl.nrm.begin(), afl.nrm.end());
		tau.insert(tau.begin(), afl.tau.begin(), afl.tau.end());
		len.insert(len.begin(), afl.len.begin(), afl.len.end());
		viscousStress = afl.viscousStress;
		lowLeft = afl.lowLeft;
		upRight = afl.upRight;
		v = afl.v;
		gammaThrough.insert(gammaThrough.begin(), afl.gammaThrough.begin(), afl.gammaThrough.end());
	};
	
	/// Деструктор
	virtual ~Airfoil() { };
	
	/// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
	void CalcNrmTauLen();	

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
	void Rotate(double alpha);			 
	
	/// \brief Масштабирование профиля 
	///
	/// Масштабирует профиль в factor раз относительно центра масс
	/// \param[in] factor масштабный коэффициент
	void Scale(double factor);			

	/// \brief Перемещение профиля 
	///
	/// Плоскопараллельно перемещает профиль на вектор \f$ \overrightarrow{dr} \f$
	///
	/// \param[in] dr константная ссылка на вектор перемещения
	void Move(const Point2D& dr);
	
	/// \brief Вычисляет габаритный прямоугольник профиля
	///
	/// Заполняет значения полей lowLeft и upRight по габаритному прямоугольнику профиля с учетом зазора
	///
	/// \param[in] gap относительная величина зазора в долях от размера габаритного прямоугольника (по умолчанию 0.02, что соответствует 2 %)
	void GetGabarits(double gap = 0.02); 
	
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
};

#endif

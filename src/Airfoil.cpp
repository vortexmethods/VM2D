/*!
\file
\brief Файл кода с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Airfoil.h"
#include "StreamParser.h"


// Конструктор
Airfoil::Airfoil()
	:np(0)
{ }


// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
void Airfoil::CalcNrmTauLen()		
{
	if (nrm.size() != np)
	{
		nrm.resize(np);
		tau.resize(np);
		len.resize(np);
	}

	Point2D rpan;

	for (size_t i = 0; i < np; ++i)
	{
		rpan = (r[i + 1] - r[i]);
		len[i] = rpan.length();
		tau[i] = rpan.unit();
		nrm[i] = { tau[i][1], -tau[i][0] };
	}
}//CalcNrmTauLen()


//Проверка, идет ли вершина i следом за вершиной j
bool Airfoil::isAfter(int i, int j) const
{
	return ((i == j + 1) || (i == 0 && j == np - 1));
}//isAfter(...)


//Перемещение профиля 
void Airfoil::Move(const Point2D& dr)	//перемещение профиля как единого целого на вектор dr
{
	for (size_t i = 0; i < np + 1; i++)
		r[i] += dr;
	GetGabarits();
}//Move(...)


//Поворот профиля 
void Airfoil::Rotate(double alpha)	//поворот профиля на угол alpha вокруг центра масс
{
	numvector<numvector<double, 2>, 2> rotMatrix = { { cos(alpha), -sin(alpha) }, { sin(alpha), cos(alpha) } };

	for (size_t i = 0; i < np + 1; i++)
		r[i] = rcm + dot(rotMatrix, r[i] - rcm);

	CalcNrmTauLen();
	GetGabarits();
}//Rotate(...)


//Масштабирование профиля
void Airfoil::Scale(double factor)	//масштабирование профиля на коэффициент factor относительно центра масс
{
	for (size_t i = 0; i < np + 1; i++)
		r[i] = rcm + factor*(r[i] - rcm);

	CalcNrmTauLen(); //строго говоря, меняются при этом только длины, нормали и касательные - без изменения
	GetGabarits();
}//Scale(...)


//Вычисляет габаритный прямоугольник профиля
void Airfoil::GetGabarits(double gap)	//определение габаритного прямоугольника
{
	lowLeft = {  1E+10,  1E+10 };
	upRight = { -1E+10, -1E+10 };

	for (size_t i = 0; i < np + 1; i++)
	{
		lowLeft[0] = std::min(lowLeft[0], r[i][0]);
		lowLeft[1] = std::min(lowLeft[1], r[i][1]);

		upRight[0] = std::max(upRight[0], r[i][0]);
		upRight[1] = std::max(upRight[1], r[i][1]);
	}

	Point2D size = upRight - lowLeft;
	lowLeft -= size*gap;
	upRight += size*gap;
}//GetGabarits(...)


//Определяет, находится ли точка с радиус-вектором r внутри габаритного прямоугольника профиля
bool Airfoil::isInsideGabarits(const Point2D& r) const
{
	return (r[0] <= upRight[0] && (r[0] >= lowLeft[0] && r[1] >= lowLeft[1] && r[1] <= upRight[1]));
}//isInsideGabarits(...)


//Определяет, находится ли точка с радиус-вектором r вне габаритного прямоугольника профиля
bool Airfoil::isOutsideGabarits(const Point2D& r) const
{
	return (r[0] > upRight[0] || (r[0] < lowLeft[0] || r[1] < lowLeft[1] || r[1] > upRight[1]));
}//isOutsideGabarits(...)

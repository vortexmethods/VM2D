/*!
\file
\brief В РАБОТЕ Заголовочный файл с описанием класса VelocityFourier
\warning Пока в состоянии нерабочего прототипа
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VELOCITYFOURIER_H
#define VELOCITYFOURIER_H

#include "unsupported/Eigen/FFT"

#include "Velocity.h"

/*!
\brief В РАБОТЕ Класс, определяющий способ вычисления скоростей

\warning Пока в состоянии нерабочего прототипа

Способ вычисления скоростей
- при помощи решения вспомогательного уравнения для функции тока на грубой сетке, используя быстрое преобразование Фурье

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/
class VelocityFourier : public Velocity
{
private:
	numvector<int, 2> nNode;
	Point2D L;
	Point2D corner;

	numvector<int, 2> nCell;
	Point2D h;

	double(*Kernel)(double);

	std::vector<std::vector<int>> vic;  //списки вихрей, попавших в ячейки сетки
	
	int cellIj(int i, int j) const; //пересчет глобального номера ячейки по "локальным координатам" ячейки
	
	//double W(const Point2D& pt) const; //"размазанная" завихренность

	void FillVic();  //заполнение списков в массиве vic

	double w(int i, int j) const; //завихренность, приведенная к узлу (i,j)

	Point2D Skos(const Point2D& R, const Point2D& X) const
	{
		double dst2 = std::max(dist2(R, X), 1e-10);
		Point2D res = { -(R[1] - X[1]), (R[0] - X[0]) };
		res *= IDPI / dst2;
		return res;
	}

public:
	VelocityFourier(const Parallel& parallel_, const Wake& wake_) :
		Velocity(parallel_, wake_)
	{
		nNode = { 16, 16 };
		L = { 1.0, 1.0 };
		corner = { 0.0, 0.0 };

		Kernel = M4;

		nCell = { nNode[0] - 1, nNode[1] - 1 };
		h = { L[0] / nCell[0], L[1] / nCell[1] };
	};

	virtual ~VelocityFourier() {};

	virtual void CalcConvVelo(double dt = 1.0);

	
};




#endif
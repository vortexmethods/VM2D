/*!
\file
\brief Файл кода с описанием класса Point2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "Point2D.h"


MPI_Datatype Point2D::mpiPoint2D;


//Конструктор и приведение типа из numvector<double, 2>
Point2D::Point2D(const numvector<double, 2>& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
}//Point2D(...)


//Конструктор копирования
Point2D::Point2D(const Point2D& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
}//Point2D(...)


//Конструктор инициализации списком
Point2D::Point2D(const std::initializer_list<double>& z)
{
	for (size_t i = 0; i < 2; ++i)
		r[i] = *(z.begin() + i);
}//Point2D(...)


//Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
Point2D Point2D::rotated(const double angle) const
{
	Point2D res;
	double cosa = cos(angle);
	double sina = sin(angle);

	res[0] = r[0] * cosa - r[1] * sina;
	res[1] = r[0] * sina + r[1] * cosa;
	return res;
}//rotated(...)


// Cоздание MPI-описателя типа
void Point2D::CreateMpiType()
{
	int          len[3] = { 1, 2, 1 };
	MPI_Aint     pos[3] = { 0, offsetof(Point2D, r), sizeof(Point2D) };
	MPI_Datatype typ[3] = { MPI_LB, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(3, len, pos, typ, &mpiPoint2D);
	MPI_Type_commit(&mpiPoint2D);
}//CreateMpiType()

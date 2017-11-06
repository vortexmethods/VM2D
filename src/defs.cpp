/*!
\file
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "defs.h"


//Сохранение матрицы в поток
void SaveToStream(const Eigen::MatrixXd& matr, std::ostream& str)
{
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение комплекснозначной матрицы в поток
void SaveToStream(const Eigen::MatrixXcd& matr, std::ostream& str)
{
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение вектора в поток
void SaveToStream(const Eigen::VectorXd& vec, std::ostream& str)
{
	for (int i = 0; i < vec.size(); ++i)
		str << vec(i) << " ";
}//SaveToStream(...)


//Сохранение списка из двумерных векторов (точек) в поток
void SaveToStream(const std::vector<Point2D>& vec, std::ostream& str)
{
	for (size_t i = 0; i < vec.size(); ++i)
		str << "{ " << vec[i][0] << " " << vec[i][1] << " } ";
}//SaveToStream(...)


//Ядро сглаживания (Монагана)
double M4(double t)
{
	double t2 = t*t;
	double mt = t > 0 ? t : -t;

	return (mt > 2.0) ? 0.0 : \
		(mt > 1.0) ? 0.5*sqr(2.0 - mt)*(1.0 - mt) : 1.0 - 2.5*t2 + 1.5*t2*mt;
}//M4(...)
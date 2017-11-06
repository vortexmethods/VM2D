/*!
\file
\brief Файл кода с описанием класса BoundaryConstLayerAver
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "BoundaryConstLayerAver.h"


//Возврат размерности вектора решения 
int BoundaryConstLayerAver::GetUnknownsSize() const
{
	return afl->np;
}//GetUnknownsSize()


//Пересчет решения на интенсивность вихревого слоя
void BoundaryConstLayerAver::SolutionToFreeVortexSheet(const Eigen::VectorXd& sol)
{
	for (size_t j = 0; j < afl->np; ++j)
		sheets.freeVortexSheet[j][0] = sol(j);
}//SolutionToFreeVortexSheet(...)


//Генерация блока матрицы
void BoundaryConstLayerAver::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl->np;

	//Panel vectors
	std::vector<Point2D> dd;
	for (size_t i = 0; i < np; ++i)
		dd.push_back(afl->tau[i] * afl->len[i]);

	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = afl->len[i];
	}

	//auxillary scalars
	numvector<double, 3> alpha, lambda;

	//auxillary vectors
	Point2D p1, s1, p2, s2, i00;
	numvector<Point2D, 3> v;

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		if (i != j)
		{
			const Point2D& di = dd[i];
			const Point2D& dj = dd[j];

			const Point2D& taui = afl->tau[i];
			const Point2D& tauj = afl->tau[j];

			p1 = CC[i + 1] - CC[j + 1];
			s1 = CC[i + 1] - CC[j];
			p2 = CC[i] - CC[j + 1];
			s2 = CC[i] - CC[j];

			alpha = { \
				afl->isAfter(j, i) ? 0.0 : Alpha(s2, s1), \
				Alpha(s2, p1), \
				afl->isAfter(i, j) ? 0.0 : Alpha(p1, p2) \
			};

			lambda = { \
				afl->isAfter(j, i) ? 0.0 : Lambda(s2, s1), \
				Lambda(s2, p1), \
				afl->isAfter(i, j) ? 0.0 : Lambda(p1, p2) \
			};

			v = { Omega(s1, taui, tauj), -Omega(di, taui, tauj), Omega(p2, taui, tauj) };

			i00 = IDPI / afl->len[i] * ((alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2]) + (lambda[0] * v[0] + lambda[1] * v[1] + lambda[2] * v[2]).kcross());

			matr(i, j) = i00 * afl->tau[i];
		}
	}

	for (size_t i = 0; i < np; ++i)
	{
		matr(i, i) = -0.5;
	}	
}//FillMatrixSelf(...)


//Заполнение правой части
void BoundaryConstLayerAver::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs)
{
	size_t np = afl->np;
	for (size_t i = 0; i < np; ++i)
		rhs(i) = -(V0*afl->tau[i]);
}//FillRhs(...)

/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.5    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/02/20     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: nummatrix.h                                                      |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Описание класса nummatrix
\author Марчевский Илья Константинович
\version 1.5   
\date 20 февраля 2019 г.
*/

#ifndef NUMMATRIX_H_
#define NUMMATRIX_H_

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Шаблонный класс, определяющий матрицу фиксированного размера
	\n Фактически представляет собой массив, для которого определено большое количество различных операций.
	\n Для доступа к элементам матрицы используется оператор [][]
	\n Оператор[] имитирует возврат ссылки/константной ссылки на numvector

	\tparam T тип элементов матрицы
	\tparam n число строк
	\tparam m число столбцов

	\author Марчевский Илья Константинович
	\version 1.5
	\date 20 февраля 2019 г.
	*/


	template<typename T, size_t n, size_t m>
	class nummatrix
	{
	public:

		/// Собственно, вектор с данными
		T r[n*m];
	public:

		/// \brief Перегрузка оператора "[]" доступа к строке
		///
		/// \tparam T тип данных
		/// \param[in] i номер строки, к которой происходит обращение
		/// \return Ссылка на numvector, соответствующий i-й строке матрицы
		numvector<T, m>& operator[](size_t i)
		{
			return *(reinterpret_cast<numvector<T, m>*>(r + (i*m)));
		}

		/// \brief Перегрузка оператора "[]" доступа к строке
		///
		/// \tparam T тип данных
		/// \param[in] i номер строки, к которой происходит обращение
		/// \return Константная сылка на numvector, соответствующий i-й строке матрицы
		const numvector<T, m>& operator[](size_t i) const
		{
			return *(reinterpret_cast<const numvector<T, m>*>(r + (i*m)));
		}


		// \brief Перегрузка оператора "=" присваивания
		//
		// \tparam T тип данных
		// \tparam n длина вектора
		// \param[in] vec константная ссылка на присваиваемую матрицу
		// \return ссылка на саму себя
		//template<typename P>
		//nummatrix<T, n, m>& operator=(const nummatrix<P, n, m>& vec)
		//{
		//	for (size_t i = 0; i < n; ++i)
		//	for (size_t j = 0; j < m; ++j)
		//		r[i*m+j] = vec[i][j];
		//	return *this;
		//}//operator=


		/// \brief Перегрузка оператора "*=" домножения матрицы на действительное число 
		///
		/// \tparam T тип данных
		/// \tparam P тип множителя
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] c числовой множитель типа, приводимого к типу компонент матрицы
		/// \return ссылка на саму себя после домножения на число
		template <typename P>
		nummatrix<T, n, m>& operator*=(const P c)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] *= c;
			return *this;
		}//operator*=


		/// \brief Перегрузка оператора "/=" деления матрицы на действительное число
		///
		/// \tparam T тип данных
		/// \tparam P тип множителя
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] c числовой делитель типа, приводимого к типу компонент матрицы
		/// \return ссылка на саму себя после деления на число
		template <typename P>
		nummatrix<T, n, m>& operator/=(const P c)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] /= c;
			return *this;
		}//operator/=


		/// \brief Перегрузка оператора "+=" прибавления другой матрицы
		///
		/// \tparam T тип данных
		/// \tparam P тип данных прибавляемой матрицы
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на прибавляемую матрицу
		/// \return ссылка на саму себя после сложения с другой матрицей
		template <typename P>
		nummatrix<T, n, m>& operator+=(const nummatrix<P, n, m>& y)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] += y.r[i];
			return *this;
		}//operator+=


		/// \brief Перегрузка оператора "-=" вычитания другой матрицы
		///
		/// \tparam T тип данных
		/// \tparam P тип данных вычитаемой матрицы
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на вычитаемую матрицу
		/// \return ссылка на саму себя после вычитания другой матрицы
		template <typename P>
		nummatrix<T, n, m>& operator-=(const nummatrix<P, n, m>& y)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] -= y.r[i];
			return *this;
		}//operator-=


		/// \brief Перегрузка оператора "+" сложения двух матриц
		///
		/// \tparam T тип данных
		/// \tparam P тип данных матрицы - второго слагаемого
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на прибавляемую матрицу
		/// \return результат сложения двух матриц, приведенный к нужному типу
		template <typename P>
		auto operator+(const nummatrix<P, n, m>& y) const -> nummatrix<typename std::remove_const<decltype(r[0] + y.r[0])>::type, n, m>
		{
			nummatrix<typename std::remove_const<decltype(r[0] + y.r[0])>::type, n, m> res;
			for (size_t i = 0; i < n*m; ++i)
				res[0][i] = r[i] + y.r[i];
			return res;
		}//operator+


		/// \brief Перегрузка оператора "-" вычитания двух матриц
		///
		/// \tparam T тип данных
		/// \tparam P тип данных матрицы - вычитаемого
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на вычитаемую матрицу
		/// \return результат вычитания двух матриц, приведенный к нужному типу
		template <typename P>
		auto operator-(const nummatrix<P, n, m>& y) const -> nummatrix<typename std::remove_const<decltype(r[0] - y.r[0])>::type, n, m>
		{
			nummatrix<typename std::remove_const<decltype(r[0] - y.r[0])>::type, n, m> res;
			for (size_t i = 0; i < n*m; ++i)
				res[0][i] = r[i] - y.r[i];
			return res;
		}//operator-


		/// \brief Перегрузка оператора "*" умножения матрицы справа на число
		///
		/// \tparam T тип данных
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] с число-множитель
		/// \return результат умножения вектора на число, приведенный к соответствующему типу
		template <typename P>
		auto operator*(const P c) const -> nummatrix<typename std::remove_const<decltype(r[0] * c)>::type, n, m>
		{
			nummatrix<typename std::remove_const<decltype(r[0] * c)>::type, n, m> res;
			for (size_t i = 0; i < n*m; ++i)
				res[0][i] = c * r[i];
			return res;
		}//operator*



		/// \brief Перегрузка оператора "-" унарного минуса
		///
		/// \tparam T тип данных
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \return противоположная матрица
		nummatrix<T, n, m> operator-() const
		{
			nummatrix<T, n, m> res;
			for (size_t i = 0; i < n*m; ++i)
				res[0][i] = -r[i];
			return res;
		}//operator-


		/// \brief Перегрузка оператора "==" логического равенства
		///
		/// \tparam T тип данных
		/// \tparam P тип данных матрицы, с которой производится сравнение
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на сравниваемую матрицу
		/// \return true, если матрицы одинаковые, false в противном случае
		template <typename P>
		bool operator==(const nummatrix<P, n, m>& y) const
		{
			for (size_t i = 0; i < n*m; ++i)
				if (r[i] != y.r[i])
					return false;
			return true;
		}//operator==


		/// \brief Перегрузка оператора "!=" логического неравенства
		///
		/// \tparam T тип данных
		/// \tparam P тип данных матрицы, с которой производится сравнение
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] y константная ссылка на сравниваемую матрицу
		/// \return true, если матрицы различаются, false в противном случае
		template <typename P>
		bool operator!=(const nummatrix<P, n, m>& y) const
		{
			return !(*this == y);
		}//operator!=


		/// \brief Вычисление размерности матрицы (числа строк и столбцов в ней)
		///
		/// \return размерность вектора (число элементов в нем)
		std::pair<size_t, size_t> size() const { return{ n, m }; }


		/// \brief Вычисление 1-нормы матрицы
		///
		/// \return 1-норма матрицы
		auto norm1() const -> typename std::remove_const<typename std::remove_reference<decltype(r[0])>::type>::type
		{
			std::vector<typename std::remove_const<std::remove_reference<decltype(r[0])>::type>::type> ms(m, 0);
			for (size_t i = 0; i < m; ++i)
				for (size_t j = 0; j < n; ++j)
					ms[i] += abs(r[j*m + i]);
			return std::max_element(ms.begin(), ms.end());
		}

		/// \brief Вычисление inf-нормы матрицы
		///
		/// \return inf-норма матрицы
		auto norminf() const -> typename std::remove_const<typename std::remove_reference<decltype(r[0])>::type>::type
		{
			std::vector<typename std::remove_const<typename std::remove_reference<decltype(r[0])>::type>::type> ms(n, 0);
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j)
					ms[i] += abs(r[i*m + j]);
			return std::max_element(ms.begin(), ms.end());
		}


		/// Установка всех компонент матрицы в константу (по умолчанию --- нуль)
		/// \tparam T тип данных
		/// \tparam n число строк
		/// \tparam m число столбцов
		///
		/// \param[in] val константа, значению которой приравниваются все компоненты матрицы (по умолчанию 0.0)
		/// \return ссылка на саму матрицу	
		nummatrix<T, n, m>& toZero(T val = 0.0)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] = val;
			return *this;
		}


		/// Установка матрицы в единичную 
		/// \tparam T тип данных
		/// \tparam n число строк и столбцов
		///
		/// \return ссылка на саму матрицу	
		nummatrix<T, n, n>& toIdentity()
		{
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j)
					r[i*m + j] = (i == j) ? 1 : 0;
			return *this;
		}


		/// Пустой конструктор
		nummatrix() { };


		/// \brief Конструктор, инициализирующий всю матрицу одной и той же константой
		///
		/// \tparam P тип данных инициализирующей константы
		/// \param[in] z значение, которым инициализируются все компоненты матрицы
		template <typename P>
		explicit nummatrix(const T z)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] = z;
		}//nummatrix(...)


		/// \brief Конструктор копирования
		///
		/// \tparam P тип данных копируемой матрицы
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] mtr константная ссылка на копируемую матрицу
		template <typename P>
		nummatrix(const nummatrix<P, n, m>& mtr)
		{
			for (size_t i = 0; i < n*m; ++i)
				r[i] = mtr.r[i];
		}//nummatrix(...)


		/// \brief Конструктор инициализации с помощью std::vector из std::vector
		///
		/// \tparam P тип данных инициализирующего std::vector
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] vec константная ссылка на инициализирующий вектор
		template <typename P>
		nummatrix(const std::vector<std::vector<P>>& vec)
		{
			if ((vec.size() != n) || (vec[0].size() != m))
				throw;
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j)
					r[i*m + j] = vec[i][j];
		}//nummatrix(...)


		/// \brief Конструктор инициализации при помощи numvector из numvector
		///
		/// \tparam P тип данных в инициализирующем numvector
		/// \tparam n число строк
		/// \tparam m число столбцов
		/// \param[in] mtr константная ссылка на numvector из numvector
		template <typename P>
		nummatrix(const numvector<numvector<P, m>, n>& mtr)
		{
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j)
					r[i*m + j] = mtr[i][j];
		}//nummatrix(...)


#if !defined(__CUDACC__)
	/// \brief Нешаблонный конструктор инициализации списком из numvector
	///
	/// \tparam T тип данных
	/// \param[in] z константная ссылка на список инициализации
	/// \warning Длина списка инициализации проверяется, при несовпадении бросается исключение
		nummatrix(const std::initializer_list<numvector<T, m>>& z)
		{
			if (z.size() != n)
				throw;
			for (size_t i = 0; i < n; ++i)
			{
				const numvector<T, m>& numv = *(z.begin() + i);
				for (size_t j = 0; j < m; ++j)
					r[i*m + j] = numv[j];
			}
		}//nummatrix(...)


		/// \brief Шаблонный конструктор инициализации списком из numvector
		///
		/// \tparam P тип данных векторов в списке инициализации
		/// \param[in] z константная ссылка на список инициализации
		/// \warning Длина списка инициализации проверяется, при несовпадении бросается исключение
		template <typename P>
		nummatrix(const std::initializer_list<numvector<P, m>>& z)
		{
			if (z.size() != n)
				throw;
			for (size_t i = 0; i < n; ++i)
			{
				const numvector<T, m>& numv = *(z.begin() + i);
				for (size_t j = 0; j < m; ++j)
					r[i*m + j] = numv[j];
			}
		}//nummatrix(...)
#endif 


	/// \brief Приведение матрицы к типу std::vector из std::vector
	///
	/// \return вектор типа std::vector, состоящий из std::vector, являющихся строками матрицы
		operator std::vector<std::vector<T>>() const
		{
			std::vector<std::vector<T>> vec;
			vec.reserve(n);
			for (size_t i = 0; i < n; ++i)
			{
				vec[i].reserve(m);
				for (size_t j = 0; j < m; ++j)
					vec[i].push_back(r[i*m + j]);
			}
			return vec;
		}


		/// Транспонирование
		inline nummatrix<T, m, n> transpose() const
		{
			nummatrix<T, m, n> res;
			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < m; ++j)
					res[j][i] = r[i*m + j];
			return res;
		}//Transpose(...)


		/// \brief Умножение матрицы на вектор
		///
		/// \tparam P тип данных вектора
		/// \tparam n число строк матрицы
		/// \tparam m число столбцов матрицы и длина вектора 
		/// \param[in] x константная ссылка на вектор
		/// \return вектор --- результат умножения матрицы на вектор
		template<typename P>
		auto operator&(const numvector<P, m>& x) -> \
			numvector<typename std::remove_const<decltype(r[0] * x[0])>::type, n>
		{
			numvector<typename std::remove_const<decltype(r[0] * x[0])>::type, n> res;

			for (size_t i = 0; i < n; ++i)
				res[i] = (*this)[i] & x;

			return res;
		}//operator&

		/// \brief Умножение матрицы на матрицу
		///
		/// \tparam P тип данных матрицы - второго множителя
		/// \tparam n число строк матрицы
		/// \tparam m число столбцов матрицы и длина вектора 
		/// \param[in] B константная ссылка на матрицу	
		/// \return вектор --- результат умножения матрицы на вектор
		template<typename P, size_t p>
		auto operator&(const nummatrix<P, m, p>& B) -> \
			nummatrix<typename std::remove_const<decltype(r[0] * x[0])>::type, n, p>
		{
			nummatrix<typename std::remove_const<decltype(r[0] * x[0])>::type, n, p> res(0);

			for (size_t i = 0; i < n; ++i)
				for (size_t j = 0; j < n; ++j)
					for (size_t k = 0; k < n; ++k)
						res[i][j] += (*this)[i][k] * B[k][j];
			return res;
		}//operator&

	};


	/// \brief Умножение вектора на матрицу
	///
	/// \tparam T тип данных матрицы
	/// \tparam P тип данных вектора
	/// \tparam n длина вектора и число строк матрицы
	/// \tparam m число столбцов матрицы
	/// \param[in] A константная ссылка на матрицу
	/// \param[in] x константная ссылка на вектор
	/// \return вектор --- результат умножения матрицы на вектор
	template<typename T, typename P, size_t n, size_t m>
	auto operator&(const numvector<P, n>& x, const nummatrix<T, n, m>& A) -> \
		numvector<typename std::remove_const<decltype(x[0] * A.r[0])>::type, m>
	{
		numvector<typename std::remove_const<decltype(x[0] * A.r[0])>::type, m> res;
		for (size_t i = 0; i < m; ++i)
		{
			res[i] = 0;
			for (size_t j = 0; j < n; ++j)
				res[i] += x[j] * A[j][i];
		}
		return res;
	}//operator&(...)


	/// \brief Умножение вектора на вектор внешним образом
	template<typename T, typename P, size_t n, size_t m>
	inline auto operator| (const numvector<T, n>& x, const numvector<P, m>& y) -> \
		nummatrix<typename std::remove_const<decltype(x[0] * y[0])>::type, n, m>
	{
		nummatrix<typename std::remove_const<decltype(x[0] * y[0])>::type, n, m> res;

		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < m; ++j)
				res[i][j] = x[i] * y[j];

		return res;
	}//operator|(...)



	/// \brief Перегрузка оператора "<<" вывода в поток
	///
	/// Выводит в поток вектор, элементы которого записаны в фигурных скобках и разделены запятой с пробелом 
	///
	/// \tparam T тип данных
	/// \tparam n размерность вектора
	/// \param[in,out] str ссылка на поток вывода
	/// \param[in] x константная ссылка на вектор
	/// \return ссылка на поток вывода
	template<typename T, size_t n, size_t m>
	std::ostream& operator<< (std::ostream& str, const nummatrix<T, n, m>& x)
	{
		str << "{ ";
		for (size_t j = 0; j < n - 1; ++j)
			str << x[j] << ", ";
		str << x[n - 1];
		str << " }";
		return str;
	}//operator<<

}//namespace VMlib

using VMlib::nummatrix;

#endif
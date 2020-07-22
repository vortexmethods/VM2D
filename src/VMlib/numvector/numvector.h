/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.9    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/07/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: numvector.h                                                      |
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
\brief Описание класса numvector
\author Марчевский Илья Константинович
\version 1.9   
\date 22 июля 2020 г.
*/

#ifndef NUMVECTOR_H_
#define NUMVECTOR_H_

/*
#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : warning: " #desc)
*/

/*
#if defined(__GNUC__) || defined(__clang__)
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif
*/

#define DEPRECATED

#include <array>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <ostream>
#include <set>
#include <vector>

namespace VMlib
{

	class Point2D;
	class Point;

	/*!
	\brief Шаблонный класс, определяющий вектор фиксированной длины
	\n Фактически представляет собой массив, для которого определено большое количество различных операций.
	\n Для доступа к элементам массива используется оператор []

	\tparam T тип элементов вектора
	\tparam n длина вектора

	\author Марчевский Илья Константинович
	\version 1.9
	\date 22 июля 2020 г.
	*/
		

	template<typename T, size_t n>
	class numvector : public std::array<T, n>
	{
	protected:

	public:


		/// \brief Оператор "&" скалярного умножения
		///
		/// \tparam P тип данных компонент вектора - второго сомножителя
		/// \tparam n длина обоих векторов
		/// \param[in] y константная ссылка на второй множитель
		/// \return результат вычисления скалярного произведения, приведенный к нужному типу	
		template <typename P>
		auto operator& (const numvector<P, n>& y) const -> typename std::remove_const<decltype(this->data()[0] * y[0])>::type
		{
			typename std::remove_const<decltype(this->data()[0] * y[0])>::type res = 0;
			for (size_t j = 0; j < n; ++j)
				res += this->data()[j] * y[j];
			return res;
		}//operator&(...)


		/// \brief Оператор "^" векторного произведения 
		///
		/// Определен только для трехмерных векторов
		///
		/// \tparam P тип данных компонент вектора - второго сомножителя
		/// \param[in] y константная ссылка на второй множитель
		/// \return результат вычисления векторного произведения, приведенный к нужному типу
		template <typename P>
		auto operator^(const numvector<P, 3>& y) const -> numvector<typename std::remove_const<decltype(this->data()[1] * y[2])>::type, 3>
		{
			numvector<typename std::remove_const<decltype(this->data()[1] * y[2])>::type, 3> vec;
			vec[0] = this->data()[1] * y[2] - this->data()[2] * y[1];
			vec[1] = this->data()[2] * y[0] - this->data()[0] * y[2];
			vec[2] = this->data()[0] * y[1] - this->data()[1] * y[0];
			return vec;
		}//operator^(...)


		/// \brief Оператор "^" вычисления третьей компоненты векторного произведения
		///
		/// Определен только для двумерных векторов
		///
		/// \tparam P тип данных компонент вектора - второго множителя
		/// \param[in] y константная ссылка на второй множитель
		/// \return результат вычисления третьей компоненты векторного произведения двух двумерных векторов, приведенный к нужному типу
		template <typename P>
		auto operator^ (const numvector<P, 2>& y) const -> typename std::remove_const<decltype(this->data()[0] * y[1])>::type
		{
			return (this->data()[0] * y[1] - this->data()[1] * y[0]);
		}//operator^(...)


		/// \brief Оператор "*=" домножения вектора на действительное число 
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam P тип данных множителя
		/// \tparam n длина вектора
		/// \param[in] c числовой множитель типа, приводимого к типу компонент вектора
		/// \return ссылку на самого себя после домножения на число
		template <typename P>
		numvector<T, n>& operator*=(P c)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] *= c;
			return *this;
		}//operator*=(...)


		/// \brief Оператор "/=" деления вектора на действительное число
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam P тип данных множителя
		/// \tparam n длина вектора
		/// \param[in] c числовой делитель типа, приводимого к типу компонент вектора
		/// \return ссылку на самого себя после деления на число
		template <typename P>
		numvector<T, n>& operator/=(P c)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] /= c;
			return *this;
		}//operator/=(...)


		/// \brief Оператор "+=" прибавления другого вектора
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam P тип данных компонент прибавляемого вектора
		/// \tparam n длина вектора
		/// \param[in] y константная ссылка на прибавляемый вектор
		/// \return ссылку на самого себя после сложения с другим вектором
		template <typename P>
		numvector<T, n>& operator+=(const numvector<P, n>& y)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] += y[i];
			return *this;
		}//operator+=(...)


		/// \brief Оператор "-=" вычитания другого вектора
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam P тип данных компонент вычитаемого вектора
		/// \tparam n длина вектора
		/// \param[in] y константная ссылка на вычитаемый вектор
		/// \return ссылка на самого себя после вычитания другого вектора
		template <typename P>
		numvector<T, n>& operator-=(const numvector<P, n>& y)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] -= y[i];
			return *this;
		}//operator-=(...)


		/// \brief Оператор "+" сложения двух векторов
		///
		/// \tparam P тип данных компонент вектора - второго слагаемого
		/// \tparam n длина вектора
		/// \param[in] y константная ссылка на прибавляемый вектор
		/// \return результат сложения двух векторов, приведенный к нужному типу
		template <typename P>
		auto operator+(const numvector<P, n>& y) const -> numvector<typename std::remove_const<decltype(this->data()[0] + y[0])>::type, n>
		{
			numvector<typename std::remove_const<decltype(this->data()[0] + y[0])>::type, n> res;
			for (size_t i = 0; i < n; ++i)
				res[i] = this->data()[i] + y[i];
			return res;
		}//operator+(...)


		/// \brief Оператор "-" вычитания двух векторов
		///
		/// \tparam P тип данных компонент вектора - вычитаемого
		/// \tparam n длина вектора
		/// \param[in] y константная ссылка на вычитаемый вектор
		/// \return результат вычитания двух векторов, приведенный к нужному типу
		template <typename P>
		auto operator-(const numvector<P, n>& y) const -> numvector<typename std::remove_const<decltype(this->data()[0] - y[0])>::type, n>
		{
			numvector<typename std::remove_const<decltype(this->data()[0] - y[0])>::type, n> res;
			for (size_t i = 0; i < n; ++i)
				res[i] = this->data()[i] - y[i];
			return res;
		}//operator-(...)


		/// \brief Оператор "*" умножения вектора на число (вектор слева, число справа)
		///
		/// \tparam P тип данных множителя
		/// \tparam n длина вектора
		/// \param[in] c число-множитель
		/// \return результат умножения вектора на число, приведенный к соответствующему типу
		template <typename P>
		auto operator*(const P c) const -> numvector<typename std::remove_const<decltype(this->data()[0] * c)>::type, n>
		{
			numvector<typename std::remove_const<decltype(this->data()[0] * c)>::type, n> res;
			for (size_t i = 0; i < n; ++i)
				res[i] = c * this->data()[i];
			return res;
		}//operator*(...)


		/// \brief Оператор "-" унарного минуса
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam n длина вектора	
		/// \return противоположный вектор
		numvector<T, n> operator-() const
		{
			numvector<T, n> res;
			for (size_t i = 0; i < n; ++i)
				res[i] = -this->data()[i];
			return res;
		}//operator-()


		/// \brief Оператор "+" унарного плюса
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam n длина вектора	
		/// \return константную ссылку на самого себя
		const numvector<T, n>& operator+() const
		{
			return *this;
		}//operator+()


		/// \brief Оператор "==" логического равенства
		///
		/// \tparam P тип данных компонент вектора, с которым производится сравнение
		/// \tparam n длина вектора	
		/// \param[in] y константная ссылка на сравниваемый вектор
		/// \return true, если векторы одинаковые, false в противном случае
		template <typename P>
		bool operator==(const numvector<P, n>& y) const
		{
			for (size_t i = 0; i < n; ++i)
				if (this->data()[i] != y[i])
					return false;
			return true;
		}//operator==(...)


		/// \brief Перегрузка оператора "!=" логического неравенства
		///
		/// \tparam P тип данных компонент вектора, с которым производится сравнение
		/// \tparam n длина вектора	
		/// \param[in] y константная ссылка на сравниваемый вектор
		/// \return true, если векторы различаются, false в противном случае
		template <typename P>
		bool operator!=(const numvector<P, n>& y) const
		{
			return !(*this == y);
		}//operator!=(...)




		/// \brief Вычисление 1-нормы вектора
		///
		/// Сумма модулей компонент вектора
		///
		/// \return 1-норма вектора
		auto norm1() const -> typename std::remove_const<typename std::remove_reference<decltype(this->data()[0])>::type>::type
		{
			typename std::remove_const<typename std::remove_reference<decltype(this->data()[0])>::type>::type res = 0;
			for (size_t i = 0; i < n; ++i)
				res += abs(this->data()[i]);
			return res;
		}//norm1()


		/// \brief Вычисление inf-нормы вектора
		///
		/// Наибольшая по модулю компонента вектора
		///
		/// \return inf-норма вектора
		auto norminf() const -> typename std::remove_const<typename std::remove_reference<decltype(this->data()[0])>::type>::type
		{
			typename std::remove_const<typename std::remove_reference<decltype(this->data()[0])>::type>::type res = 0;
			for (size_t i = 0; i < n; ++i)
			{
				if (abs(this->data()[i]) > res)
					res = abs(this->data()[i]);
			}
			return res;
		}//norminf()


		/// \brief Вычисление 2-нормы (длины) вектора
		///
		/// Корень из скалярного квадрата вектора
		///
		/// \tparam P тип результата (по умолчанию double)
		///
		/// \return норма (длина) вектора
		template <typename P = double>
		P length() const
		{
			P res = *this & *this;
			return sqrt(res);
		}//length()


		/// \brief Вычисление квадрата нормы (длины) вектора
		///
		/// Скалярный квадрат вектора
		///
		/// \return квадрат нормы (длины) вектора того же типа, что и компоненты вектора
		auto length2() const -> typename std::remove_const<typename std::remove_reference<decltype(this->data()[0])>::type>::type
		{
			return (*this & *this);
		}//length2()


		/// \brief Вычисление орта вектора или вектора заданной длины, коллинеарного данному
		///
		/// Если в качестве новой длины указано отрицательное число --- вектор будет противоположно направленным
		///
		/// \tparam P тип числа, задающего длину вектора
		/// \param[in] newlen длина получаемого вектора (по умолчанию 1.0)
		/// \return вектор, коллинеарный исходному, заданной длины (по умолчанию 1.0)
		///
		/// \warning для получения float-орта от вектора с компонентами типа float или целыми нужно явно указать параметр 1.0f
		template <typename P = double>
		auto unit(P newlen = 1) const -> numvector<typename std::remove_const<decltype(this->data()[0] * newlen)>::type, n>
		{
			auto ilen = static_cast<decltype(this->data()[0] * newlen)>(newlen / std::max(this->length(), 1e-16));
			return (*this * ilen);
		}//unit(...)



		/// \brief Нормирование вектора на заданную длину
		///
		/// Если в качестве новой длины указано отрицательное число --- у вектора будет изменено направление
		///
		/// \tparam P тип числа, задающего длину вектора
		/// \param[in] newlen новая длина вектора (по умолчанию 1.0)
		///
		/// \warning Работает корректно только для векторов с компонентами типа float и double
		template <typename P = double>
		void normalize(P newlen = 1.0)
		{
			auto ilen = static_cast<decltype(this->data()[0] * newlen)>(newlen / std::max(this->length(), 1e-16));
			*this *= ilen;
		}//normalize(...)



		/// \brief Проверка вхождения элемента в вектор
		///
		/// \tparam P тип данных проверяемого элемента
		/// \param[in] s проверяемый элемент
		/// \return позиция первого вхождения элемента s; если не входит --- возвращает (-1), приведенный к типу size_t
		template <typename P>
		size_t member(const P& s) const
		{
			for (size_t i = 0; i < n; ++i)
				if (this->data()[i] == s)
					return i;

			return static_cast<size_t>(-1);
		}//member(...)



		/// \brief Приведение вектора к типу std::set
		///
		/// \tparam P тип данных компонент множества
		/// \return множество типа std::set, состоящее из тех же элементов, что исходный вектор
		template <typename P>
		operator std::set<P>() const
		{
			std::set<P> newset;
			for (size_t i = 0; i < n; ++i)
				newset.insert(this->data()[i]);
			return newset;
		}//toSet()


		/// \brief Приведение вектора к типу std::vector
		///
		/// \tparam P тип данных компонент std::vector
		/// \return вектор типа std::vector, состоящий из тех же элементов, что исходный вектор
		template <typename P>
		operator std::vector<P>() const
		{
			std::vector<P> vec;
			vec.reserve(n);
			for (size_t i = 0; i < n; ++i)
				vec.push_back(this->data()[i]);
			return vec;
		}


		/// \brief "Вращение" вектора на несколько позиций влево
		///
		/// Исходный вектор при этом не изменяется
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam n длина вектора	
		/// \param[in] k количество позиций, на которые производится "вращение"
		/// \return вектор, полученный "вращением" исходного на k позиций влево
		numvector<T, n> rotateLeft(size_t k) const
		{
			if (k >= n)
				throw;
			numvector<T, n> res;

			//for (size_t i = 0; i < n; ++i) 
			//	res[i] = r[(i + k) % n];

			for (size_t i = 0; i < n - k; ++i)
				res[i] = this->data()[i + k];
			for (size_t i = n - k + 1; i < n; ++i)
				res[i] = this->data()[i + k - n];
			return res;
		}//rotateLeft(...)



		/// \brief Геометрический поворот двумерного вектора на 90 градусов
		///
		/// Исходный вектор при этом не изменяется
		/// \n Эквивалентно умножению слева на орт третьей оси, т.е. \f$ \vec k \times \vec r \f$
		///
		/// \tparam T тип данных
		/// \return новый двумерный вектор, полученный поворотом исходного на 90 градусов
		numvector<T, 2> kcross() const
		{
			numvector<T, 2> res;
			res[0] = -this->data()[1];
			res[1] =  this->data()[0];
			return res;
		}//kcross()


		/// \brief Установка всех компонент вектора в константу (по умолчанию --- нуль)
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam P тип данных константы
		/// \tparam n длина вектора
		/// \param[in] val константа, значению которой приравниваются все компоненты вектора (по умолчанию 0)
		/// \return ссылка на сам вектор
		template <typename P = T>
		numvector<T, n>& toZero(P val = 0)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = val;
			return *this;
		}


		/// \brief Вычисление квадрата расстояния до другой точки
		///
		/// \tparam P тип данных второй точки	
		/// \param[in] y константная ссылка на радиус-вектор второй точки
		/// \return квадрат расстояния между точками
		template<typename P>
		auto dist2To(const numvector<P, n>& y) const -> typename std::remove_const<decltype(this->data()[0] - y[0])>::type
		{
			return (*this - y).length2();
		}//dist2To(...)


		/// \brief Вычисление расстояния между двумя точками
		///
		/// \tparam P тип данных второй точки
		/// \tparam R тип данных результата (по умолчанию double)
		/// \param[in] y константная ссылка на радиус-вектор второй точки
		/// \return расстояние между точками
		template<typename R = double, typename P>
		R distTo(const numvector<P, n>& y)
		{
			R res = (*this - y) & (*this - y);
			return sqrt(res);
		}//distTo(...)


		/// Пустой конструктор
		numvector() { };


		/// \brief Конструктор, инициализирующий весь вектор одной и той же константой
		///
		/// \tparam P тип данных присваиваемой константы
		/// \param[in] c значение, которым инициализируются все компоненты вектора
		template <typename P>
		explicit numvector(const P c)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = c;
		}//numvector(...)
		


		/// \brief Нешаблонный конструктор копирования
		///
		/// \tparam T тип данных компонент вектора
		/// \tparam n длина вектора
		/// \param[in] vec константная ссылка на копируемый вектор
		numvector(const numvector<T, n>& vec)
		{
			//for (size_t i = 0; i < n; ++i)
			//	r[i] = vec[i];
			memcpy(this->data(), vec.data(), n * sizeof(T));
		}//numvector(...)


		/// \brief Шаблонный конструктор копирования
		///
		/// \tparam P тип данных копируемого вектора
		/// \tparam n длина вектора
		/// \param[in] vec константная ссылка на копируемый вектор
		template <typename P>
		numvector(const numvector<P, n>& vec)
		{
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = vec[i];
		}//numvector(...)


		/// \brief Конструктор инициализации с помощью std::vector
		///
		/// \tparam P тип данных инициализирующего std::vector
		/// \tparam n длина вектора
		/// \param[in] vec константная ссылка на инициализирующий вектор
		///
		/// \warning при несовпадении длины кидает исключение
		template <typename P>
		numvector(const std::vector<P>& vec)
		{
			if (vec.size() != n)
				throw;
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = vec[i];
		}//numvector(...)


#if !defined(__CUDACC__)
	/// \brief Нешаблонный конструктор инициализации списком
	///
	/// \tparam T тип данных
	/// \param[in] z константная ссылка на список инициализации
	///
	/// \warning Длина списка инициализации проверяется, при несовпадении бросается исключение
		numvector(const std::initializer_list<T>& z)
		{
			if (z.size() != n)
				throw;
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = *(z.begin() + i);
		}//numvector(...)



		/// \brief Шаблонный конструктор инициализации списком
		///
		/// \tparam P тип данных инициализирующего списка
		/// \param[in] z константная ссылка на список инициализации
		///
		/// \warning Длина списка инициализации проверяется, при несовпадении бросается исключение
		template <typename P>
		numvector(const std::initializer_list<P>& z)
		{
			if (z.size() != n)
				throw;
			for (size_t i = 0; i < n; ++i)
				this->data()[i] = *(z.begin() + i);
		}//numvector(...)
#endif 


	/// \brief Явный конструктор инициализации вектором другой размерности
	///
	/// - если новая размерность меньше старой --- лишние элементы "отбрасываются"
	/// - если новая размерность больше старой --- новые элементы заполняются заданным элементом
	///
	/// \tparam P тип данных инициализирующего вектора
	/// \tparam p размерность инициализирующего вектора
	/// \param[in] vec константная ссылка на инициализирующий вектор
	/// \param[in] add элемент для заполнения при необходимости недостающих компонент (по умолчанию 0)
		template <typename P, size_t p>
		explicit numvector(const numvector<P, p>& vec, T add = 0)
		{
			size_t minPN = (p < n) ? p : n;
			for (size_t i = 0; i < minPN; ++i)
				this->data()[i] = vec[i];
			for (size_t i = minPN; i < n; ++i)
				this->data()[i] = add;
		}//numvector(...)	


//		////////////////////////////////////////////////////////////////////////////
//		//// Далее deprecate-функции для ПОЛНОЙ СОВМЕСТИМОСТИ со старой версией ////
//		////                             --------------------                   ////
//		////////////////////////////////////////////////////////////////////////////
//
//		/// \brief Оператор присваивания всем компонентам вектора одного и того же числа
//		DEPRECATED numvector<T, n>& operator=(double c)
//		{
//			toZero(c);
//			return *this;			
//		}//operator(...)


//		/// \brief Построение множества std::set на основе вектора
//		///
//		/// \return множество типа std::set, состоящее из тех же элементов, что исходный вектор
//		DEPRECATED std::set<T> toSet() const
//		{
//			//(deprecate: use implicit type conversion)
//
//			std::set<T> newset;
//			for (size_t i = 0; i < n; ++i)
//				newset.insert(this->data()[i]);
//			return newset;
//		}//toSet()


//		/// \brief Построение вектора std::vector на основе вектора
//		///
//		/// \return вектор типа std::vector, состоящий из тех же элементов, что исходный вектор
//		DEPRECATED std::vector<T> toVector() const
//		{
//			//(deprecate: use implicit type conversion)
//
//			std::vector<T> vec;
//			vec.reserve(n);
//			for (size_t i = 0; i < n; ++i)
//				vec.push_back(this->data()[i]);
//			return vec;
//		}

	}; //class numvector

//    /// \todo Исследовать целесообразность наличия нешаблонного умножения
//	inline numvector<double, 3> operator*(double c, const numvector<double, 3>& x)
//	{
//		numvector<double, 3> res(x);
//		for (size_t i = 0; i < 3; ++i)
//			res[i] *= c;
//		return res;
//	}//operator*(...)
	
	/// \brief Оператор "*" умножения вектора на число (число слева, вектор справа)
	///
	/// \tparam T тип данных компонент вектора и множителя	
	/// \tparam n длина вектора
	/// \param[in] c числовой множитель
	/// \param[in] x константная ссылка на умножаемый вектор
	/// \return результат умножения вектора на число, приведенный к соответствующему типу
	template<typename T, size_t n>
	inline numvector<T, n> operator*(double c, const numvector<T, n>& x)
	{
		numvector<T, n> res(x);
		for (size_t i = 0; i < n; ++i)
			res[i] *= c;
		return res;
	}//operator*(...)


	/// \brief Оператор "*" умножения вектора на число (число слева, вектор справа)
	///
	/// \tparam T тип данных компонент вектора
	/// \tparam P тип данных числового множителя
	/// \tparam n длина вектора
	/// \param[in] c числовой множитель
	/// \param[in] x константная ссылка на умножаемый вектор
	/// \return результат умножения вектора на число, приведенный к соответствующему типу
	template<typename T, typename P, size_t n>
	auto operator*(const P c, const numvector<T, n>& x) -> numvector<typename std::remove_const<decltype(x[0] * c)>::type, n>
	{
		numvector<typename std::remove_const<decltype(x[0] * c)>::type, n> res;
		for (size_t i = 0; i < n; ++i)
			res[i] = x[i] * c;
		return res;
	}//operator*(...)


#if !defined(__CUDACC__)
 /// \brief Быстрое вычисление векторного произведения
 ///
 /// Определено только для трехмерных векторов
 /// \n Оптимизировано за счет отсутствия вызова конструктора, предполагает наличие трех уже созданных векторов
 ///
 /// \tparam T тип данных компонент вектора первого множителя
 /// \tparam P тип данных компонент вектора второго множителя
 /// \tparam R тип данных компонент вектора результата
 /// \param[in] x константная ссылка на первый множитель
 /// \param[in] y константная ссылка на второй множитель
 /// \param[out] z ссылка на результат векторного умножения
	template<typename T, typename P, typename R>
	inline void cross(const numvector<T, 3>& x, const numvector<P, 3>& y, numvector<R, 3>& z)
	{
		z = { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
	}//cross(...)
#endif


/// \brief Вычисление квадрата расстояния между двумя точками
///
/// \tparam T тип данных компонент первого радиус-вектора
/// \tparam P тип данных компонент второго радиус-вектора
/// \tparam n размерность векторов
/// \param[in] x константная ссылка на радиус-вектор первой точки
/// \param[in] y константная ссылка на радиус-вектор второй точки
/// \return квадрат расстояния между точками
	template<typename T, typename P, size_t n>
	inline auto dist2(const numvector<T, n>& x, const numvector<P, n>& y) -> typename std::remove_const<decltype(x[0] - y[0])>::type
	{
		numvector<typename std::remove_const<decltype(x[0] - y[0])>::type, n> p = x - y;
		return p.length2();
	}//dist2(...)


	/// \brief Вычисление расстояния между двумя точками
	///
	/// \tparam T тип данных компонент радиус-вектора первой точки
	/// \tparam P тип данных компонент радиус-вектора второй точки
	/// \tparam R тип данных результата (по умолчанию double)
	/// \param[in] x константная ссылка на радиус-вектор первой точки
	/// \param[in] y константная ссылка на радиус-вектор второй точки
	/// \return расстояние между точками
	template<typename R = double, typename T, typename P, size_t n>
	inline R dist(const numvector<T, n>& x, const numvector<P, n>& y)
	{
		numvector<R, n> p = x - y;
		return sqrt(p&p);
	}//dist(...)


	/// \brief Оператор "<<" вывода вектора в поток
	///
	/// Выводит в поток вектор, элементы которого записываются в фигурных скобках и разделены запятой с пробелом 
	///
	/// \tparam T тип данных компонент вектора
	/// \tparam n размерность вектора
	/// \param[in,out] str ссылка на поток вывода
	/// \param[in] x константная ссылка на вектор
	/// \return ссылку на поток вывода
	template<typename T, size_t n>
	std::ostream& operator<< (std::ostream& str, const numvector<T, n>& x)
	{
		str << "{ ";
		for (size_t j = 0; j < n - 1; ++j)
			str << x[j] << ", ";
		str << x[n - 1];
		str << " }";
		return str;
	}//operator<<(...)


	////////////////////////////////////////////////
	//// Далее --- операциями с парами векторов ////
	////////////////////////////////////////////////


	/// \brief Оператор прибавления "+=" для пар векторов
	/// 
	/// tparam T тип данных вектора - первого компонента в паре
	/// tparam P тип данных вектора - второго компонента в паре
	/// tparam R тип данных вектора - первого компонента в прибавляемой паре
	/// tparam S тип данных вектора - второго компонента в прибавляемой паре
	/// tparam n длина векторов - компонентов обеих пар
	/// param[in] a ссылка на первую пару
	/// param[in] b константная ссылка на прибавляемую пару
	/// return ссылку на первую пару после прибавления к ней второй пары
	template<typename T, typename P, typename R, typename S, size_t n>
	inline std::pair<numvector<T, n>, numvector<P, n>>& operator+=(std::pair<numvector<T, n>, numvector<P, n>>& a, const std::pair<numvector<R, n>, numvector<S, n>>& b)
	{
		a.first += b.first;
		a.second += b.second;
		return a;
	}//operator+=(...)


	/// \brief Оператор сложения "+" для пар векторов
	/// 
	/// tparam T тип данных вектора - первого компонента в левой паре
	/// tparam P тип данных вектора - второго компонента в левой паре
	/// tparam R тип данных вектора - первого компонента в правой паре
	/// tparam S тип данных вектора - второго компонента в правой паре
	/// tparam n длина векторов - компонентов обеих пар
	/// param[in] a константная ссылка на левую пару
	/// param[in] b константная ссылка на правую пару
	/// return пару, получаемую при покомпонентном суммировании левой и правой пар, приведенную к нужному типу
	template<typename T, typename P, typename R, typename S, size_t n>
	inline auto operator+(const std::pair<numvector<T, n>, numvector<P, n>>& a, const std::pair<numvector<R, n>, numvector<S, n>>& b) -> \
		std::pair<numvector<typename std::remove_const<decltype(a.first[0] + b.first[0])>::type, n>, numvector<typename std::remove_const<decltype(a.second[0] + b.second[0])>::type, n >>
	{
		std::pair<numvector<typename std::remove_const<decltype(a.first[0] + b.first[0])>::type, n>, numvector<typename std::remove_const<decltype(a.second[0] + b.second[0])>::type, n >> res;
		res.first = a.first + b.first;
		res.second = a.second + b.second;
		return res;
	}//operator+(...)


	/// \brief Оператор домножения "*=" пары векторов на число (пара слева, число справа)
	/// 
	/// tparam T тип данных вектора - первого компонента в паре
	/// tparam P тип данных вектора - второго компонента в паре
	/// tparam R тип данных числового множителя
	/// tparam n длина векторов - компонентов пары
	/// param[in] a ссылка на пару
	/// param[in] c числовой множитель
	/// return ссылка на пару после ее покомпонентного домножения на число
	template<typename T, typename P, typename R, size_t n>
	inline std::pair<numvector<T, n>, numvector<P, n>>& operator*=(std::pair<numvector<T, n>, numvector<P, n>>& a, R c)
	{
		a.first *= c;
		a.second *= c;
		return a;
	}//operator*=(...)


	/// \brief Оператор умножения "*" числа на пару векторов (число слева, пара справа)
	/// 
	/// tparam T тип данных вектора - первого компонента в паре
	/// tparam P тип данных вектора - второго компонента в паре
	/// tparam R тип данных числового множителя
	/// tparam n длина векторов - компонентов пары
	/// param[in] c числовой множитель
	/// param[in] a константная ссылка на пару
	/// return пару, получаемую при покомпонентном умножении пары на число, приведенную к нужному типу
	template<typename T, typename P, typename R, size_t n>
	inline auto operator*(R c, const std::pair<numvector<T, n>, numvector<P, n>>& a) -> \
		std::pair<numvector<typename std::remove_const<decltype(c * a.first[0])>::type, n>, numvector<typename std::remove_const<decltype(c * a.second[0])>::type, n >>
	{
		std::pair<numvector<typename std::remove_const<decltype(c * a.first[0])>::type, n>, numvector<typename std::remove_const<decltype(c * a.second[0])>::type, n >> res;
		res.first = c * a.first;
		res.second = c * a.second;
		return res;
	}//operator*(...)


	/// \brief Оператор умножения "*" пары векторов на число (пара слева, число справа)
	/// 
	/// tparam T тип данных вектора - первого компонента в паре
	/// tparam P тип данных вектора - второго компонента в паре
	/// tparam R тип данных числового множителя
	/// tparam n длина векторов - компонентов пары
	/// param[in] c числовой множитель
	/// param[in] a константная ссылка на пару
	/// return пару, получаемую при покомпонентном умножении пары на число, приведенную к нужному типу
	template<typename T, typename P, typename R, size_t n>
	inline auto operator*(const std::pair<numvector<T, n>, numvector<P, n>>& a, R c) -> \
		std::pair<numvector<typename std::remove_const<decltype(a.first[0] * c)>::type, n>, numvector<typename std::remove_const<decltype(a.second[0] * c)>::type, n >>
	{
		return c * a;
	}//operator*(...)


	/// \brief Оператор "<<" вывода пары векторов в поток
	///
	/// Выводит в поток векторы из пары, разделяя их запятой с пробелом 
	///
	/// \tparam T тип данных вектора - первого компонента пары
	/// \tparam P тип данных вектора - второго компонента пары
	/// \tparam n размерность векторов
	/// \param[in,out] str ссылка на поток вывода
	/// \param[in] x константная ссылка на вектор
	/// \return ссылку на поток вывода
	template<typename T, typename P, size_t n>
	std::ostream& operator<< (std::ostream& str, const std::pair<numvector<T, n>, numvector<P, n>>& x)
	{
		str << "{ " << x.first << ", " << x.second << " }";
		return str;
	}//operator<<


	////////////////////////////////////////////////////////////////////////
	// Далее deprecate-функции для ПОЛНОЙ СОВМЕСТИМОСТИ со старой версией //
	//                             --------------------                   //
	////////////////////////////////////////////////////////////////////////

//	/// \brief Умножение квадратной матрицы на вектор (без приведения типов)
//	///
//	/// \tparam T тип данных компонент матрицы и вектора
//	/// \tparam n размерность матрицы и вектора
//	/// \param[in] A константная ссылка на матрицу
//	/// \param[in] x константная ссылка на вектор
//	/// \return вектор результат умножения матрицы на вектор
//	template<typename T, size_t n>
//	DEPRECATED inline numvector<T, n> dot(const numvector<numvector<T, n>, n>& A, const numvector<T, n>& x)
//	{
//		//deprecate: use nummatrix.operator& instead of dot(numvector<numvector<>>, numvector<>)
//
//		numvector<T, n> res;
//		for (size_t i = 0; i < n; ++i)
//			res[i] = A[i] & x;
//		return res;
//	}//dot(...)


//	/// \brief Умножение вектора на квадратную матрицу (без приведения типов)
//	///
//	/// \tparam T тип данных компонент матрицы и вектора
//	/// \tparam n размерность матрицы и вектора
//	/// \param[in] A константная ссылка на матрицу
//	/// \param[in] x константная ссылка на вектор
//	/// \return вектор результат умножения матрицы на вектор
//	template<typename T, size_t n>
//	DEPRECATED inline numvector<T, n> dot(const numvector<T, n>& x, const numvector<numvector<T, n>, n>& A)
//	{
//		//deprecate: use operator&(numvector<>, nummatrix<>) instead of dot(numvector<>, numvector<numvector<>>)
//
//		numvector<T, n> res;
//		for (size_t i = 0; i < n; ++i)
//		{
//			res[i] = 0.0;
//			for (size_t j = 0; j < n; ++j)
//				res[i] += x[j] * A[j][i];
//		}
//		return res;
//	}//dot(...)


//	/// \brief Умножение квадратной матрицы на вектор
//	///
//	/// \tparam T тип данных вектора и матрицы
//	/// \tparam n длина вектора и размерность матрицы
//	/// \param[in] A константная ссылка на матрицу
//	/// \param[in] x константная ссылка на вектор
//	/// \return вектор --- результат умножения матрицы на вектор
//	template<typename T, size_t n>
//	DEPRECATED inline numvector<T, n> matDotVec(const numvector<numvector<T, n>, n>& A, const numvector<T, n>& x)
//	{
//		//deprecate: use nummatrix.operator& instead of matDotVec(numvector<numvector<>>, numvector<>)
//
//		numvector<T, n> res;
//		for (size_t i = 0; i < n; ++i)
//			res[i] = A[i] * x;
//		return res;
//	}//dot(...)


//	/// \brief Умножение вектора на вектор той же размерности внешним образом
//	///
//	/// \tparam T тип данных компонент векторов
//	/// \tparam n размерность векторов
//	/// \param[in] x константная ссылка на первый вектор
//	/// \param[in] y константная ссылка на второй вектор
//	/// \return матрицу ранга 1, являющюуся внешним (кронекеровым) произведением двух векторов
//	template<typename T, size_t n>
//	DEPRECATED inline numvector<numvector<T, n>, n> KronProd(const numvector<T, n>& x, const numvector<T, n>& y)
//	{
//		//deprecate: use nummatrix.operator|(numvector<>, numvector<>) instead of KronProd(numvector<>, numvector<>)
//		numvector<numvector<T, n>, n> res;
//		for (size_t i = 0; i < n; ++i)
//			for (size_t j = 0; j < n; ++j)
//				res[i][j] = x[i] * y[j];
//		return res;
//	}//KronProd(...)


//	/// \brief Транспонирование квадратной матрицы
//	///
//	/// \tparam T тип данных компонент матрицы
//	/// \tparam n размерность матрицы
//	/// \return транспонированную матрицу
//	template<typename T, size_t n>
//	DEPRECATED inline numvector<numvector<T, n>, n> Transpose(const numvector<numvector<T, n>, n>& A)
//	{
//		//deprecate: use member class nummatrix.transpose() instead of Transpose(numvector<numvector>)
//
//		numvector<numvector<T, n>, n> res;
//		for (size_t i = 0; i < n; ++i)
//			for (size_t j = 0; j < n; ++j)
//				res[i][j] = A[j][i];
//		return res;
//	}


	/// \brief Вычисление третьей компоненты векторного произведения
	///
	/// Определено для векторов любой размерности, используются только первые 2 компоненты исходных векторов
	///
	/// \tparam T тип данных компонент вектора - первого множителя
	/// \tparam P тип данных компонент вектора - второго множителя
	/// \tparam n размерность исходных векторов
	/// \param[in] x константная ссылка на первый множитель
	/// \param[in] y константная ссылка на второй множитель
	/// \return третья компонента вектороного произведения
	template<typename T, typename P, size_t n>
	DEPRECATED inline double cross3(const numvector<T, n>& x, const numvector<P, n>& y)
	{
		//deprecate: use numvector.operator^ instead of cross3(numvector<>,numvector<>)

		return (x[0] * y[1] - x[1] * y[0]);
	}//cross3(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const Point2D& x, const Point2D& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> xx = *(reinterpret_cast<const numvector<double, 2>*>(&x));
//		const numvector<double, 2> yy = *(reinterpret_cast<const numvector<double, 2>*>(&y));
//		return (xx & yy);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const numvector<double, 2>& x, const Point2D& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> yy = *(reinterpret_cast<const numvector<double, 2>*>(&y));
//		return (x & yy);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const Point2D& x, const numvector<double, 2>& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> xx = *(reinterpret_cast<const numvector<double, 2>*>(&x));
//		return (xx & y);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const Point& x, const Point& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> xx = *(reinterpret_cast<const numvector<double, 2>*>(&x));
//		const numvector<double, 2> yy = *(reinterpret_cast<const numvector<double, 2>*>(&y));
//		return (xx & yy);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const numvector<double, 2>& x, const Point& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> yy = *(reinterpret_cast<const numvector<double, 2>*>(&y));
//		return (x & yy);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const Point& x, const numvector<double, 2>& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		const numvector<double, 2> xx = *(reinterpret_cast<const numvector<double, 2>*>(&x));
//		return (xx & y);
//	}//operator*(...)




//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const numvector<double, 2>& x, const numvector<double, 2>& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		return (x & y);
//	}//operator*(...)


//	/// \brief Скалярное умножение двух векторов
//	/// 
//	/// \param[in] x константная ссылка на первый множитель
//	/// \param[in] y константная ссылка на второй множитель
//	/// \return результат скалярного умножения
//	DEPRECATED inline double operator*(const numvector<double, 3>& x, const numvector<double, 3>& y)
//	{
//		//deprecate: use numvector.operator& instead of operator*
//
//		return (x & y);
//	}//operator*(...)

}//namespace VMlib

using VMlib::numvector;
using VMlib::dist2;

#endif
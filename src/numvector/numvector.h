/*!
\file
\brief Описание класса numvector
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef NUMVECTOR_H_
#define NUMVECTOR_H_

#include <vector>
#include <set>

/*!
\brief Шаблонный класс, определяющий вектор фиксированной длины
\n Фактически представляет собой массив, для которого определено большое количество различных операций.
\n Для доступа к элементам массива используется оператор []
\tparam T тип элементов вектора
\tparam n длина вектора
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


template<typename T, int n>
class numvector
{
protected:

	/// Собственно, вектор с данными
	T r[n];
public:

	/// \brief Перегрузка оператора "[]" доступа к элементу
	///
	/// \tparam T тип данных
	/// \param[in] j номер элемента, к которому происходит обращение
	/// \return Ссылка на элемент
	T& operator[](int j) { return r[j]; }


	/// \brief Перегрузка оператора "[]" доступа к элементу
	///
	/// \tparam T тип данных
	/// \param[in] j номер элемента, к которому происходит обращение
	/// \return Константная ссылка на элемент
	const T& operator[](int j) const { return r[j]; }


	/// \brief Перегрузка оператора "=" присваивания
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] vec константная ссылка на присваиваемый вектор
	/// \return ссылка на самого себя
	numvector<T, n>& operator=(const numvector<T, n>& vec) 
	{
		for (int i = 0; i < n; ++i)
			r[i] = vec.r[i];
		return *this;
	}//operator=


	/// \brief Перегрузка оператора "*" скалярного умножения
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] y константная ссылка на второй множитель
	/// \return результат вычисления скалярного произведения
	double operator* (const numvector<T, n>& y) const 
	{
		T res = 0;
		for (int j = 0; j < n; ++j)
			res += r[j] * y[j];
		return res;
	}//operator*


	/// \brief Перегрузка оператора "^" векторного произведения 
	///
	/// Определен только для трехмерных векторов
	///
	/// \tparam T тип данных
	/// \param[in] y константная ссылка на второй множитель
	/// \return результат вычисления векторного произведения
	numvector<T, 3> operator^ (const numvector<T, 3>& y) const
	{
		return{ r[1] * y[2] - r[2] * y[1], r[2] * y[0] - r[0] * y[2], r[0] * y[1] - r[1] * y[0] };
	}//operator^
	

	/// \brief Перегрузка оператора "^" вычисления третьей компоненты векторного произведения
	///
	/// Определен только для двумерных векторов
	///
	/// \tparam T тип данных
	/// \param[in] y константная ссылка на второй множитель
	/// \return результат вычисления третьей компоненты векторного произведения двух двумерных векторов
	T operator^ (const numvector<T, 2>& y) const
	{
		return r[0] * y[1] - r[1] * y[0];
	}//operator^


	/// \brief Перегрузка оператора "*=" домножения вектора на действительное число 
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] m числовой множитель, приводится к типу double
	/// \return ссылка на самого себя после домножения на число
	numvector<T, n>& operator*=(const double m)
	{		
		for (int i = 0; i < n; ++i)
			r[i] *= m;
		return *this;
	}//operator*=


	/// \brief Перегрузка оператора "/=" деления вектора на действительное число
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] m числовой делитель, приводится к типу double
	/// \return ссылка на самого себя после деления на число
	numvector<T, n>& operator/=(const double m)
	{
		for (int i = 0; i < n; ++i)
			r[i] /= m;
		return *this;
	}//operator/=


	/// \brief Перегрузка оператора "+=" прибавления другого вектора
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] y константная ссылка на прибавляемый вектор
	/// \return ссылка на самого себя после сложения с другим вектором
	numvector<T, n>& operator+=(const numvector<T, n>& y)
	{
		for (int i = 0; i < n; ++i)
			r[i] += y[i];
		return *this;
	}//operator+=


	/// \brief Перегрузка оператора "-=" вычитания другого вектора
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] y константная ссылка на вычитаемый вектор
	/// \return ссылка на самого себя после вычитания другого вектора
	numvector<T, n>& operator-=(const numvector<T, n>& y)
	{
		for (int i = 0; i < n; ++i)
			r[i] -= y[i];
		return *this;
	}//operator-=


	/// \brief Перегрузка оператора "+" сложения двух векторов
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] y константная ссылка на прибавляемый вектор
	/// \return результат сложения двух векторов
	numvector<T, n> operator+(const numvector<T, n>& y) const
	{
		numvector<T, n> res(*this);
		for (int i = 0; i < n; ++i)
			res[i] += y[i];
		return res;
	}//operator+
	

	/// \brief Перегрузка оператора "-" вычитания двух векторов
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] y константная ссылка на вычитаемый вектор
	/// \return результат вычитания двух векторов
	numvector<T, n> operator-(const numvector<T, n>& y) const
	{
		numvector<T, n> res(*this);
		for (int i = 0; i < n; ++i)
			res[i] -= y[i];
		return res;
	}//operator-


	/// \brief Перегрузка оператора "*" умножения вектора справа на число
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] m число-множитель, приводится к типу double
	/// \return результат вычитания двух векторов
	numvector<T, n> operator*(const double m) const
	{
		numvector<T, n> res(*this);
		res *= m;
		return res;
	}//operator*


	/// \brief Перегрузка оператора "-" унарного минуса
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора	
	/// \return противоположный вектор
	numvector<T, n> operator-() const
	{
		numvector<T, n> res(0);
		for (int i = 0; i < n; ++i)
			res[i] = -(*this)[i];
		return res;
	}//operator*


	/// \brief Вычисление размерности вектора (числа элементов в нем)
	///
	/// \return размерность вектора (число элементов в нем)
	int size() const { return n; }


	/// \brief Вычисление нормы (длины) вектора
	///
	/// Корень из скалярного квадрата вектора
	///
	/// \return норма (длина) вектора
	double length() const { return sqrt(*this * *this); }

	
	/// \brief Вычисление квадрата нормы (длины) вектора
	///
	/// Скалярный квадрат вектора
	///
	/// \return квадрат нормы (длины) вектора
	double length2() const { return (*this * *this); }


	/// \brief Вычисление орта вектора или вектора заданной длины, коллинеарного данному
	///
	/// Если в качестве новой длины указано отрицательное число --- вектор будет противоположно направленным
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора	
	/// \param[in] newlen длина получаемого вектора (по умолчанию 1.0)
	///
	/// \return вектор, коллинеарный исходному, заданной длины (по умолчанию 1.0)
	numvector<T, n> unit(double newlen = 1.0) const
	{
		numvector<T, n> res(*this);
		double ilen = newlen / this->length();
		res *= ilen;
		return res;
	}//unit()


	/// \brief Нормирование вектора на заданную длину
	///
	/// Если в качестве новой длины указано отрицательное число --- у вектора будет изменено направление
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора	
	/// \param[in] newlen новая длина вектора (по умолчанию 1.0)
	///
	/// \return ссылка на самого себя после нормирования на заданную длину (по умолчанию 1.0)
	numvector<T, n>& normalize(double newlen = 1.0) 
	{
		double ilen = newlen / this->length();
		*this *= ilen;
		return (*this);
	}//normalize()


	/// \brief Проверка вхождения элемента в вектор
	///
	/// \param[in] s проверяемый элемент
	///
	/// \return позиция первого вхождения элемента s; если не входит --- возвращает (-1), приведенный к типу size_t
	size_t member(T s)
	{
		for (int i = 0; i < n; ++i)
		if (r[i] == s)
			return i;

		return (-1);
	}//member(...)


	/// \brief Построение множества std::set на основе вектора
	///
	/// \return множество типа std::set, состоящее из тех же элементов, что исходный вектор
	std::set<T> toSet() const
	{
		set<T> newset;
		for (size_t i = 0; i < n; ++i)
			newset.insert(r[i]);
		return newset;
	}//toSet()


	/// \brief "Вращение" вектора на несколько позиций влево
	///
	/// Исходный вектор при этом не изменяется
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора	
	/// \param[in] k количество позиций, на которые производится "вращение"
	///
	/// \return вектор, полученный "вращением" исходного на k позиций влево
	numvector<T, n> rotateLeft(size_t k) const
	{
		numvector<T, n> res;
		for (int i = 0; i < n; ++i)
			res[i] = r[(i + k) % n];
		return res;
	}//rotateLeft(...)

	
	/// \brief Геометрический поворот двумерного вектора на 90 градусов
	///
	/// Исходный вектор при этом не изменяется
	/// \n Эквивалентно умножению слева на орт третьей оси, т.е. \f$ \vec k \times \vec r \f$
	///
	/// \tparam T тип данных	
	///
	/// \return новый двумерный вектор, полученный поворотом исходного на 90 градусов
	numvector<T, 2> kcross() const
	{
		return { -r[1], r[0] };
	}//kcross()

	
	/// Пустой конструктор
	numvector()	{ };


	/// \brief Конструктор, инициализирующий весь вектор одной и той же константой
	///
	/// \tparam T тип данных
	/// \param[in] z значение, которым инициализируются все компоненты вектора
	numvector(const T z)
	{
		for (int i = 0; i < n; ++i)
			r[i] = z;
	}//numvector(...)
	

	/// \brief Конструктор копирования
	///
	/// \tparam T тип данных
	/// \tparam n длина вектора
	/// \param[in] vec константная ссылка на копируемый вектор
	numvector(const numvector<T, n>& vec)
	{
		for (int i = 0; i < n; ++i)
			r[i] = vec.r[i];
	}//numvector(...)


	/// \brief Конструктор инициализации списком
	///
	/// \tparam T тип данных
	/// \param[in] z константная ссылка на список инициализации
	/// \warning Длина списка инициализации не проверяется, от него берутся только необходимое число первых элементов
	numvector(const std::initializer_list<T>& z)
	{
		//if (z.size() > n) exit(0);
		//for (size_t i = 0; i < z.size(); ++i)
		//	r[i] = *(z.begin() + i); 
		//for (size_t i = z.size(); i < n; ++i)
		//	r[i] = 0;

		for (size_t i = 0; i < n; ++i)
			r[i] = *(z.begin() + i);
	}//numvector(...)

}; //class numvector


/// \brief Перегрузка оператора "*" умножения вектора слева на число
///
/// \tparam T тип данных вектора
/// \tparam n длина вектора
/// \param[in] m число-множитель, которое приводится к типу double
/// \param[in] x константная ссылка на умножаемый вектор
/// \return результат умножения вектора на число
template<typename T, int n>
inline numvector<T, n> operator*(const double m, const numvector<T, n>& x)
{
	numvector<T, n> res(x);
	res *= m;
	return res;
}//operator*


/// \brief Умножение квадратной матрицы на вектор
///
/// \tparam T тип данных вектора и матрицы
/// \tparam n длина вектора и размерность матрицы
/// \param[in] A константная ссылка на матрицу
/// \param[in] x константная ссылка на вектор
/// \return вектор --- результат умножения матрицы на вектор
template<typename T, int n>
inline numvector<T, n> dot(const numvector<numvector<T, n>, n>& A, const numvector<T, n>& x)
{
	numvector<T, n> res;

	for (int i = 0; i < n; ++i)
		res[i] = A[i] * x;
	
	return res;
}//dot(...)


/// \brief Быстрое вычисление векторного произведения
///
/// Определено только для трехмерных векторов
/// \n Оптимизировано за счет отсутствия вызова конструктора, предполагает наличие трех уже созданных векторов
///
/// \tparam T тип данных
/// \param[in] x константная ссылка на первый множитель
/// \param[in] y константная ссылка на второй множитель
/// \param[out] z ссылка на результат векторного умножения
template<typename T>
inline void cross(const numvector<T, 3>& x, const numvector<T, 3>& y, numvector<T, 3>& z) 
{
	z = { x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0] };
}//cross(...)


/// \brief Вычисление третьей компоненты векторного произведения 
///
/// Определено для векторов любой размерности, используются только первые 2 компоненты исходных векторов
///
/// \tparam T тип данных
/// \tparam n размерность исходных векторов
/// \param[in] x константная ссылка на первый множитель
/// \param[in] y константная ссылка на второй множитель
/// \return третья компонента вектороного произведения
template<typename T, int n>
inline double cross3(const numvector<T, n>& x, const numvector<T, n>& y)
{
	return (x[0] * y[1] - x[1] * y[0]);
}//cross3(...)


/// \brief Вычисление квадрата расстояния между двумя точками
///
/// \tparam T тип данных
/// \tparam n размерность векторов
/// \param[in] x константная ссылка на радиус-вектор первой точки
/// \param[in] y константная ссылка на радиус-вектор второй точки
/// \return квадрат расстояния между точками
template<typename T, int n>
inline double dist2(const numvector<T, n>& x, const numvector<T, n>& y) 
{
	T res = 0;
	for (int j = 0; j < n; ++j)
		res += (x[j] - y[j])*(x[j] - y[j]);
	return res;	
}//dist2(...)


/// \brief Вычисление расстояния между двумя точками
///
/// \tparam T тип данных
/// \tparam n размерность векторов
/// \param[in] x константная ссылка на радиус-вектор первой точки
/// \param[in] y константная ссылка на радиус-вектор второй точки
/// \return расстояние между точками
template<typename T, int n>
inline double dist(const numvector<T, n>& x, const numvector<T, n>& y) 
{
	return sqrt(dist2(x, y));	
}//dist(...)


/// \brief Перегрузка оператора "<<" вывода в поток
///
/// Выводит в поток вектор, элементы которого записаны в фигурных скобках и разделены запятой с пробелом 
///
/// \tparam T тип данных
/// \tparam n размерность вектора
/// \param[in,out] str ссылка на поток вывода
/// \param[in] x константная ссылка на вектор
/// \return ссылка на поток вывода
template<typename T, int n>
std::ostream& operator<< (std::ostream& str, const numvector<T,n>& x)
{
	str << "{ ";
	for (int j = 0; j < n - 1; ++j)
		str << x[j] << ", ";
	str << x[n - 1];
	str << " }";
	return str;
}//operator<<


/// \brief Приведение вектора к вектору другой размерности
///
/// - если новая размерность меньше старой --- лишние элементы "отбрасываются"
/// - если новая размерность больше сиарой --- новые элементы заполняются нулями
///
/// \tparam T тип данных
/// \tparam n размерность исходного вектора
/// \tparam p размерность нового вектора
/// \param[in] r константная ссылка на исходный вектор
/// \return вектор новой размерности
template<typename T, int n, int p>
numvector<T, p> toOtherLength(const numvector<T, n>& r)
{
	numvector<T, p> res;
	for (int i = 0; i < min(p, n); ++i)
		res[i] = r[i];
	for (int i = min(p, n); i < p; ++i)
		res[i] = 0;

	return res;
}//toOtherLength(...)

#endif
/*!
\file
\brief Заголовочный файл с описанием класса Sheet
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef SHEET_H
#define SHEET_H

#include <vector>


/*!
\brief Класс, опеделяющий слои на поверхности обтекаемого профиля

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Sheet
{
public:
	
	/// Список из характеристик свободного вихревого слоя на панелях
	std::vector<std::vector<double>> freeVortexSheet;

	/// Список из характеристик присоединенного вихревого слоя на панелях
	std::vector<std::vector<double>> attachedVortexSheet;

	/// Список из характеристик присоединенного слоя источников на панелях
	std::vector<std::vector<double>> attachedSourceSheet;

	/// Пустой конструктор 
	Sheet() { };

	/// Деструктор
	~Sheet() { };

	/// \brief Установка pазмерностей всех векторов и их обнуление
	///
	/// \param[in] np число панелей на профиле (внешняя размерность списков)
	/// \param[in] layerDim количество чисел, которыми характеризуются слои на каждой из панелей
	void SetLayersDim(size_t np, size_t layerDim);
};


#endif
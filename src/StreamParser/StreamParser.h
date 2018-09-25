/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: StreamParser.h                                                   |
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
\brief Заголовочный файл с описанием класса StreamParser
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#ifndef STREAMPARSER_H
#define STREAMPARSER_H

//#define VARNAME(var) #var

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "defs.h"
#include "Point2D.h"
#include "Vortex2D.h"

/*!
\brief Класс, позволяющий выполнять разбор файлов и строк с настройками и параметрами
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/
class StreamParser
{
private:
	/// Указатель на поток вывода логов
	static std::ostream* pinfo;

	/// Указатель на поток вывода ошибок
	static std::ostream* perr;

	/// \brief Данные считываемые из основного файла-паспорта расчета
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из паспорта расчета
	/// \n По умолчанию данные берутся из файла passport, но это имя может быть изменено ключом pspFile в строке параметров задачи 
	std::unordered_map<std::string, std::vector<std::string>> database;
	
	/// \brief Данные считываемые из списка умолчаний
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из файла со списком значений по умолчанию
	/// \n По умолчанию данные берутся из файла defaults
	std::unordered_map<std::string, std::vector<std::string>> defaults;

	/// \brief Данные считываемые из параметров конкретной задачи
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий данные из строки --- перечня параметров конкретной задачи в файле problems. 
	/// \n Здесь же определяются переменные, которые затем могут использоваться методом разыменовывания в файлах-паспортах конкретных задач
	std::unordered_map<std::string, std::vector<std::string>> vars;

	/// \brief Данные считываемые из перечня параметров-переключателей
	///
	/// Неупорядоченный ассоциативный контейнер, хранящий значения параметров-переключателей 
	/// \n По умолчанию данные берутся из файла switcher. 	
	std::unordered_map<std::string, std::vector<std::string>> switchers;

	/// \brief Парсинг заданного потока 
	///
	/// Разбирает поток вида key = value   или   key = { value1, value2, ..., valueN }
	/// \n Возвращает базу данных, указывается тип скобок, по умолчанию --- круглые
	/// \param[in] stream ссылка на поток, который требуется распарсить
	/// \param[out] database ссылка на заполняемую базу данных (unordered_map)
	/// \param[in] specificKey поиск только ключей из данного списка (по умолчанию --- пустой список, что означает разбор всех ключей)
	/// \param[in] replaceVars признак замены переменных из базы vars (по умолчанию false)
	/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "(" )
	/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- ")" )	
	void ParseStream(std::istream& stream, std::unordered_map<std::string, std::vector<std::string>>& database, std::vector<std::string> specificKey = {}, bool replaceVars = false, char openBracket = '(', char closeBracket = ')');

	/// \brief Замена переменных в строке их значениями
	///
	/// Находит в строке все переменные, имена которых начинаются с символа $,
	/// и заменяет их на значения, извлекаемые из базы данных vars
	/// \param[in] st строка, которая подлежит обработке
	void ReplaceVarsInString(std::string& st);


public:
	/// \brief Конструктор, принимающий четыре потока
	/// 
	/// Заполняет 4 базы данных --- database, defaults, switchers, vars
	/// \n Используется для установки основных параметров расчета
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] defaultsStream ссылка на поток по умолчанию
	/// \param[in] switchersStream ссылка на поток параметров-переключателей
	/// \param[in] varsStream ссылка на поток переменных
	/// \param[in] specificKey поиск только ключей из данного списка (по умолчанию --- пустой список, что означает разбор всех ключей)
	StreamParser(std::istream& mainStream, std::istream& defaultsStream, std::istream& switchersStream, std::istream& varsStream, std::vector<std::string> specificKey = {});
		
	/// \brief Конструктор, принимающий два потока
	/// 
	/// Заполняет 2 базы данных --- database, defaults; базы данных switchers и vars остаются пустыми
	/// \n Используется для первоначального считывания параметров задач из списка задач и формирования очереди
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] defaultsStream ссылка на поток по умолчанию
	/// \param[in] specificKey поиск только ключей из данного списка (по умолчанию --- пустой список, что означает разбор всех ключей)
	StreamParser(std::istream& mainStream, std::istream& defaultsStream, std::vector<std::string> specificKey = {});
	
	/// \brief Конструктор, принимающий один поток
	/// 
	/// Заполняет 1 базу данных --- database; базы данных defaults, switchers и vars остаются пустыми
	/// \n Используется для считывания профиля, следа, списка задач
	///
	/// \param[in] mainStream ссылка на основной поток параметров
	/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "{" )
	/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- "}" )
	StreamParser(std::istream& mainStream, char openBracket = '{', char closeBracket = '}');

	/// Деструктор
	~StreamParser() { };
	
	/// \brief Перевод строки в верхний регистр
	///
	/// \param[in] line константная ссылка на обрабатываемую строку
	/// \return строка в верхнем регистре
	static std::string UpperCase(const std::string& line);
	
	/// \brief Pазбор строки, содержащей запятые, на отдельные строки
	///
	/// Запятые внутри парных скобок игнорируются и не воспринимаются как разделители
	/// \param[in] line разбираемая строка
	/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "(" )
	/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- ")" )
	/// \return список (вектор) из отдельных строк
	static std::vector<std::string> StringToVector(std::string line, char openBracket = '(', char closeBracket = ')');
	
	/// \brief Объединение вектора (списка) из строк в одну строку
	///
	/// После каждой строки вставляется символ конца строки
	/// \param[in] _vecString константная ссылка на вектор из строк
	/// \return сплошную строку с символами конца строки между подстроками
	static std::string VectorStringToString(const std::vector<std::string>& _vecString);
	
	/// \brief Разбор строки на пару ключ-значение
	///
	/// Разбирает строку вида xxx(yyy) на пару подстрок xxx и yyy
	/// \param[in] line разбираемая строка вида xxx(ууу)
	/// \param[in] upcase признак перевода ключа в верхний регистр
	/// \return пара строк <xxx, yyy>
	static std::pair<std::string, std::string> SplitString(std::string line, bool upcase = true);

	/// \brief Установка значения параметра по умолчанию
    ///
	/// Устанавливает параметру значение по умолчанию, если оно задано в виде константной ссылки
	/// \tparam T тип параметра, значение которого устанавливается
	/// \param[in] name константная ссылка на строку --- имя устанавливаемого параметра (необходимо только для отображения в лог)
	/// \param[in] res ссылка на задаваемый параметр
	/// \param[in] defValue указатель на константу --- значение по умолчанию
	template <typename T>
	void SetDefault(const std::string& name, T& res, const T* defValue) const
	{
		if (defValue != nullptr)
		{
			res = *defValue;
			*pinfo << "parser info: parameter <" << name << " = " << res << "> set as default" << std::endl;
		}
		else
		{
			*perr << "parser error: PARAMETER " << name << " NOT FOUND " << std::endl;
			exit(-1);
		}
	}

	/// \brief Считывание вектора из двумерных точек из базы данных
	/// 
	/// Переопределение метода get() для ситывания вектора из точек (Point2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на значение, считываемое из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, std::vector<Point2D>& res, const std::vector<Point2D>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			((search_it = defaults.find(UpperCase(name))) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];

				Point2D elem;
				size_t pos = s.find(',', 1);
				std::stringstream(s.substr(1, pos)) >> elem[0];
				std::stringstream(s.substr(pos + 1, s.length() - pos)) >> elem[1];
				res.push_back(elem);
			}
		}
		//else
		//	SetDefault(name, res, defValue);
	};
		

	/// \brief Считывание вектора из вихрей из базы данных
	/// 
	/// Переопределение метода get() для ситывания вектора из вихрей (Vortex2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на вектор из вихрей, считываемый из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, std::vector<Vortex2D>& res, const std::vector<Vortex2D>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			((search_it = defaults.find(UpperCase(name))) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];

				Point2D r;
				double g;
				size_t pos1 = s.find(',', 1);
				size_t pos2 = s.find(',', pos1 + 1);
				std::stringstream(s.substr(1, pos1)) >> r[0];
				std::stringstream(s.substr(pos1 + 1, s.length() - pos1)) >> r[1];
				std::stringstream(s.substr(pos2 + 1, s.length() - pos2)) >> g;

				Vortex2D elem(r, g);

				res.push_back(elem);
			}
		}
		//else
		//	SetDefault(name, res, defValue);
	};
	
	/// \brief Считывание вектора из простых типов из базы данных
	/// 
	/// Шаблонный метод get() для считывания вектора из простых типов из базы данных
	/// \tparam T тип считываемых данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на вектор из данных, считываемый из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	template <typename T>
	void get(const std::string& name, std::vector<T>& res, const std::vector<T>* defValue = nullptr) const
	{
		res.resize(0);
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( ((search_it = database.find(UpperCase(name))) != database.end()) ||
			 ((search_it = defaults.find(UpperCase(name))) != defaults.end()))
		{
			for (size_t i = 0; i < (search_it->second).size(); ++i)
			{
				std::string s = (search_it->second)[i];

				std::stringstream ss(s);
				T elem;
				ss >> elem;
				if (typeid(elem).name() == typeid(std::string("TestString")).name())
				{
					std::string* str = reinterpret_cast<std::string*>(&elem);
					str->erase(remove(str->begin(), str->end(), '\"'), str->end());
				}
				res.push_back(elem);
			}
		}
		else
			SetDefault(name, res, defValue);
	};	
	
	/// \brief Считывание скаляра из базы данных
	/// 
	/// Шаблонный метод get() для считывания скаляра из базы данных
	/// \tparam T тип считываемых данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка данные, считываемые из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	template <typename T>
	void get(const std::string& name, T& res, const T* defValue = nullptr) const
	{
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( (( (search_it = database.find(UpperCase(name)))  != database.end()) && ((search_it->second).size() > 0)) ||
			 (( (search_it = defaults.find(UpperCase(name)))  != defaults.end()) && ((search_it->second).size() > 0))  )
		{
			if ((search_it->second).size() == 1)
			{
				std::string s = (search_it->second)[0];

				//проверка на значение-переключатель
				if (((search_it = switchers.find(UpperCase(s))) != switchers.end()) && ((search_it->second).size() > 0))
					s = search_it->second[0];
												
				std::stringstream ss(s);
				T elem;
				ss >> elem;
				if (typeid(elem).name() == typeid(std::string("TestString")).name())
				{
					std::string* str = reinterpret_cast<std::string*>(&elem);
					str->erase(remove(str->begin(), str->end(), '\"'), str->end());
				}
				res = elem;				
			}
			else
			{
				*perr << "parser error: PARAMETER " << name << " IS VECTOR " << std::endl;
				exit(-1);
			}
		}
		else
			SetDefault(name, res, defValue);
	};

	/// \brief Считывание точки из базы данных
	/// 
	/// Переопределение метода get() для ситывания точки (Point2D) из базы данных
	/// \param[in] name константная ссылка на строку --- имя считываемого параметра
	/// \param[out] res ссылка на точку, считываемую из базы данных
	/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
	void get(const std::string& name, Point2D& res, const Point2D* defValue = nullptr) const
	{
		std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it;

		//поиск ключа в базе и, если нет, то в базе умолчаний
		if ( ((search_it = database.find(UpperCase(name))) != database.end()) ||
			 ((search_it = defaults.find(UpperCase(name))) != defaults.end()) )
		{
			if ((search_it->second).size() == res.size())
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					std::stringstream ss(s);
					double elem;
					ss >> elem;
					res[i] = elem;
				}
			}
			else
			{
				*perr << "parser error: PARAMETER " << name << " LENGTH DIFFERS FROM 2" << std::endl;
				exit(-1);
			}

		}
		else
			SetDefault(name, res, defValue);
	}

};

#endif
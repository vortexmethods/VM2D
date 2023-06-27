/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.12   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2024/01/14     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2024 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: StreamParser.h                                                   |
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
\brief Заголовочный файл с описанием класса StreamParser
\author Марчевский Илья Константинович
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef STREAMPARSER_H
#define STREAMPARSER_H

//#define VARNAME(var) #var

#include <algorithm>
#include <sstream>
#include <unordered_map>

#include "defs.h"
#include "AirfoilGeometry.h"
#include "Gpu2D.h"

namespace VMlib
{
	/// \brief Функция сверки строки с шаблоном, который может содержать знаки * и ?
	///
	/// \param[in] line_ константная ссылка на проверяемую строку
	/// \param[in] pattern_ константная ссылка на строку-шаблон, с которой производится сверка
	/// rerurn признак соответствия строки шаблону
	inline int is_match(const std::string& line_, const std::string& pattern_)
	{

		const auto lineCStr = line_.c_str();
		const auto patternCStr = pattern_.c_str();

		const char* line = &(lineCStr[0]);
		const char* pattern = &(patternCStr[0]);

		// returns 1 (true) if there is a match
		// returns 0 if the pattern is not whitin the line
		int wildcard = 0;

		const char* last_pattern_start = 0;
		const char* last_line_start = 0;
		do
		{
			if (*pattern == *line)
			{
				if (wildcard == 1)
					last_line_start = line + 1;

				line++;
				pattern++;
				wildcard = 0;
			}
			else if (*pattern == '?')
			{
				if (*(line) == '\0') // the line is ended but char was expected
					return 0;
				if (wildcard == 1)
					last_line_start = line + 1;
				line++;
				pattern++;
				wildcard = 0;
			}
			else if (*pattern == '*')
			{
				if (*(pattern + 1) == '\0')				
					return 1;				

				last_pattern_start = pattern;
				//last_line_start = line + 1;
				wildcard = 1;

				pattern++;
			}
			else if (wildcard)
			{
				if (*line == *pattern)
				{
					wildcard = 0;
					line++;
					pattern++;
					last_line_start = line + 1;
				}
				else				
					line++;				
			}
			else
			{
				if ((*pattern) == '\0' && (*line) == '\0')  // end of mask
					return 1; // if the line also ends here then the pattern match
				else
				{
					if (last_pattern_start != 0) // try to restart the mask on the rest
					{
						pattern = last_pattern_start;
						line = last_line_start;
						last_line_start = 0;
					}
					else					
						return 0;					
				}
			}

		} while (*line);

		if (*pattern == '\0')
			return 1;
		else
			return 0;
	};





	/*!
	\brief Класс, позволяющий выполнять разбор файлов и строк с настройками и параметрами	
	\author Марчевский Илья Константинович
	\Version 1.12
	\date 14 января 2024 г.
	*/
	class StreamParser
	{
	private:

		/// Поток для вывода логов и сообщений об ошибках
		mutable LogStream info;

		/// \brief Данные считываемые из основного файла-паспорта расчета
		///
		/// Неупорядоченный ассоциативный контейнер, хранящий данные из паспорта расчета
		/// \n По умолчанию данные берутся из файла passport, но это имя может быть изменено ключом pspFile в строке параметров задачи 
	public:
		std::unordered_map<std::string, std::vector<std::string>> database;
	private:
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
		void ParseStream(
			std::istream& stream, 
			std::unordered_map<std::string, std::vector<std::string>>& database, 
			std::vector<std::string> specificKey = {}, 
			bool replaceVars = false, 
			char openBracket = '(', char closeBracket = ')'
		);

		/// \brief Слияние двух баз данных
		///
		/// Дополняет первую базу данных записями из второй, если такого ключа там раньше не было
		/// \n Возвращает измененную первую базу данных
		/// \param[in,out] database ссылка на дополняемую базу данных (unordered_map)
		/// \param[in] add константная ссылка на добавляемую базу данных
		void MergeDbsWithoutRepeats(std::unordered_map<std::string, std::vector<std::string>>& database, const std::unordered_map<std::string, std::vector<std::string>>&);

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
		/// \param[in, out] infoStream базовый поток для вывода логов
		/// \param[in] label константная ссылка на метку парсера для вывода логов
		/// \param[in] mainStream ссылка на основной поток параметров
		/// \param[in] defaultsStream ссылка на поток по умолчанию
		/// \param[in] switchersStream ссылка на поток параметров-переключателей
		/// \param[in] varsStream ссылка на поток переменных
		/// \param[in] specificKey поиск только ключей из данного списка (по умолчанию --- пустой список, что означает разбор всех ключей)
		StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, std::istream& defaultsStream, std::istream& switchersStream, std::istream& varsStream, std::vector<std::string> specificKey = {});

		/// \brief Конструктор, принимающий два потока
		/// 
		/// Заполняет 2 базы данных --- database, defaults; базы данных switchers и vars остаются пустыми
		/// \n Используется для первоначального считывания параметров задач из списка задач и формирования очереди
		///
		/// \param[in, out] infoStream базовый поток для вывода логов
		/// \param[in] label константная ссылка на метку парсера для вывода логов
		/// \param[in] mainStream ссылка на основной поток параметров
		/// \param[in] defaultsStream ссылка на поток по умолчанию
		/// \param[in] specificKey поиск только ключей из данного списка (по умолчанию --- пустой список, что означает разбор всех ключей)
		StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, std::istream& defaultsStream, std::vector<std::string> specificKey = {});

		/// \brief Конструктор, принимающий один поток
		/// 
		/// Заполняет 1 базу данных --- database; базы данных defaults, switchers и vars остаются пустыми
		/// \n Используется для считывания профиля, следа, списка задач
		///
		/// \param[in, out] infoStream базовый поток для вывода логов
		/// \param[in] label константная ссылка на метку парсера для вывода логов
		/// \param[in] mainStream ссылка на основной поток параметров
		/// \param[in] openBracket тип открывающейся скобки (по умолчанию --- "{" )
		/// \param[in] closeBracket тип закрывающейся скобки (по умолчанию --- "}" )
		StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, char openBracket = '{', char closeBracket = '}');

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
		/// \param[in, out] info поток для вывода логов и сообщений об ошибках
		/// \param[in] line разбираемая строка вида xxx(ууу)
		/// \param[in] upcase признак перевода ключа в верхний регистр
		/// \return пара строк (xxx, yyy)
		static std::pair<std::string, std::string> SplitString(LogStream& info, std::string line, bool upcase = true);

		/// \brief Установка значения параметра по умолчанию
		///
		/// Устанавливает параметру значение по умолчанию, если оно задано в виде константной ссылки
		/// \tparam T тип параметра, значение которого устанавливается
		/// \param[in] name константная ссылка на строку --- имя устанавливаемого параметра (необходимо только для отображения в лог)
		/// \param[in] res ссылка на задаваемый параметр
		/// \param[in] defValue указатель на константу --- значение по умолчанию
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию
		template <typename T>
		void SetDefault(const std::string& name, T& res, const T* defValue, bool echoDefault) const
		{
			if (defValue != nullptr)
			{
				res = *defValue;
				if (echoDefault)
					info('i') << "parameter <" << name << " = " << res << "> set as default" << std::endl;
			}
			else
			{
				info('e') << "parameter " << name << " is not found" << std::endl;
				exit(-1);
			}
		}

		/// \brief Считывание вектора из двумерных точек из базы данных
		/// 
		/// Переопределение метода get() для считывания вектора из точек (Point2D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на значение, считываемое из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, std::vector<Point2D>& res, const std::vector<Point2D>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;
			
			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					Point2D elem;
					size_t pos = s.find(',', 1);
					std::stringstream(s.substr(1, pos - 1)) >> elem[0];
					std::stringstream(s.substr(pos + 1, s.length() - pos - 2)) >> elem[1];
					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);
			return boolRes;
		};
		

		/// \brief Считывание вектора из трехмерных точек из базы данных
		/// 
		/// Переопределение метода get() для считывания вектора из точек (v3D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на значение, считываемое из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, std::vector<v3D>& res, const std::vector<v3D>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					v3D elem;
					size_t pos1 = s.find(',', 1);
					std::stringstream(s.substr(1, pos1 - 1)) >> elem[0];
					size_t pos2 = s.find(',', pos1 + 1);
					std::stringstream(s.substr(pos1 + 1, pos2 - pos1 - 1)) >> elem[1];
					std::stringstream(s.substr(pos2 + 1, s.length() - pos2 - 2)) >> elem[2];
					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};

		/// \brief Считывание вектора из пар чисел из базы данных
		/// 
        /// Переопределение метода get() для считывания вектора из пар чисел (std::pair) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на значение, считываемое из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		template<typename T, typename P>
		bool get(const std::string& name, std::vector<std::pair<T,P>>& res, const std::vector<std::pair<T, P>>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					std::pair<T,P> elem;
					size_t pos = s.find(',', 1);
					std::stringstream(s.substr(1, pos - 1)) >> elem.first;
					std::stringstream(s.substr(pos + 1, s.length() - pos - 2)) >> elem.second;
					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};


		/// \brief Считывание вектора из пар чисел из базы данных
		/// 
		/// Переопределение метода get() для считывания вектора из пар чисел (numvector<T,2>) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на значение, считываемое из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		template<typename T>
		bool get(const std::string& name, std::vector<numvector<T,2>>& res, const std::vector<numvector<T, 2>>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					numvector<T, 2> elem;
					size_t pos = s.find(',', 1);
					std::stringstream(s.substr(1, pos - 1)) >> elem[0];
					std::stringstream(s.substr(pos + 1, s.length() - pos - 2)) >> elem[1];
					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};


		/// \brief Считывание вектора из вихрей из базы данных
		/// 
		/// Переопределение метода get() для считывания вектора из вихрей (Vortex2D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на вектор из вихрей, считываемый из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>& res, const std::vector<Vortex2D>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					Point2D r;
					double g;
					size_t pos1 = s.find(',', 1);
					size_t pos2 = s.find(',', pos1 + 1);
					std::stringstream(s.substr(1, pos1 - 1)) >> r[0];
					std::stringstream(s.substr(pos1 + 1, pos2 - pos1 - 1)) >> r[1];
					std::stringstream(s.substr(pos2 + 1, s.length() - pos2 - 2)) >> g;

					Vortex2D elem(r, g);

					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};


		/// \brief Считывание вектора из точек, задающих профиль, из базы данных
		/// 
		/// Переопределение метода get() для считывания вектора из вихрей (Vortex2D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на вектор из вихрей, считываемый из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, std::vector<GeomPoint>& res, const std::vector<GeomPoint>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
			{
				for (size_t i = 0; i < (search_it->second).size(); ++i)
				{
					std::string s = (search_it->second)[i];

					Point2D r;
					std::string type;
					size_t pos1 = s.find(',', 1);
					size_t pos2 = s.find(',', pos1 + 1);
					std::stringstream(s.substr(1, pos1 - 1)) >> r[0];
					std::stringstream(s.substr(pos1 + 1, pos2 - pos1 - 1)) >> r[1];
					std::stringstream(s.substr(pos2 + 1, s.length() - pos2 - 2)) >> type;
					
					GeomPoint elem(r, type);


					res.push_back(elem);
				}
				boolRes = true;
			}
			//else
			//	SetDefault(name, res, defValue, echoDefault);
			

			return boolRes;
		};




		/// \brief Считывание вектора из простых типов из базы данных
		/// 
		/// Шаблонный метод get() для считывания вектора из простых типов из базы данных
		/// \tparam T тип считываемых данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на вектор из данных, считываемый из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		template <typename T>
		bool get(const std::string& name, std::vector<T>& res, const std::vector<T>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			res.resize(0);
			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
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
				boolRes = true;
			}
			else
				SetDefault(name, res, defValue, echoDefault);
			
			return boolRes;
		};


		/// \brief Считывание пары (строка, скаляр) из базы данных
		/// 
		/// Шаблонный метод get() для считывания пары (строка, скаляр) из базы данных
		/// \tparam T тип считываемых данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка данные, считываемые из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		template <typename T>
		bool get(const std::string& name, std::pair<std::string, T>& res, const std::pair<std::string, T>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if ((((search_it = database.find(UpperCase(name))) != database.end()) && ((search_it->second).size() > 0)) ||
			//	(((search_it = defaults.find(UpperCase(name))) != defaults.end()) && ((search_it->second).size() > 0)))
			if ((((search_it = search_ita) != database.end()) && ((search_it->second).size() > 0)) ||
				(((search_it = search_itb) != defaults.end()) && ((search_it->second).size() > 0)))
			{
				if ((search_it->second).size() == 1)
				{
					std::string s = (search_it->second)[0];
					res.first = s;

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
					res.second = elem;
					boolRes = true;
				}
				else
				{
					info('e') << "parameter " << name << " is list (only scalar is available)" << std::endl;
					exit(-1);
				}
			}
			else
				SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};




		/// \brief Считывание скаляра из базы данных
		/// 
		/// Шаблонный метод get() для считывания скаляра из базы данных
		/// \tparam T тип считываемых данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка данные, считываемые из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		template <typename T>
		bool get(const std::string& name, T& res, const T* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if ((((search_it = database.find(UpperCase(name))) != database.end()) && ((search_it->second).size() > 0)) ||
			//	(((search_it = defaults.find(UpperCase(name))) != defaults.end()) && ((search_it->second).size() > 0)))
			if ((((search_it = search_ita) != database.end()) && ((search_it->second).size() > 0)) ||
				(((search_it = search_itb) != defaults.end()) && ((search_it->second).size() > 0)))
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
					boolRes = true;
				}
				else
				{
					info('e') << "parameter " << name << " is list (only scalar is available)" << std::endl;
					exit(-1);
				}
			}
			else
				SetDefault(name, res, defValue, echoDefault);
			
			return boolRes;
		};

			   

		/// \brief Считывание точки из базы данных
		/// 
		/// Переопределение метода get() для считывания точки (Point2D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на точку, считываемую из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, Point2D& res, const Point2D* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
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
					boolRes = true;
				}
				else
				{
					info('e') << "parameter " << name << " length differs from 2 (only Point2D is available)" << std::endl;
					exit(-1);
				}

			}
			else
				SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		}


		/// \brief Считывание точки из базы данных
		/// 
		/// Переопределение метода get() для считывания точки (v3D) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на точку, считываемую из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, v3D& res, const v3D* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if (((search_it = database.find(UpperCase(name))) != database.end()) ||
			//	((search_it = defaults.find(UpperCase(name))) != defaults.end()))
			if (((search_it = search_ita) != database.end()) ||
				((search_it = search_itb) != defaults.end()))
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
					boolRes = true;
				}
				else
				{
					info('e') << "parameter " << name << " length differs from 3 (only v3D is available)" << std::endl;
					exit(-1);
				}

			}
			else
				SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		}


		/// \brief Считывание логической переменной из базы данных
		/// 
		/// Переопределение метода get() для ситывания логической переменной (bool) из базы данных
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка на логическую переменную, считываемую из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, bool& res, const bool* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;

			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;

			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if ((((search_it = database.find(UpperCase(name))) != database.end()) && ((search_it->second).size() > 0)) ||
			//	(((search_it = defaults.find(UpperCase(name))) != defaults.end()) && ((search_it->second).size() > 0)))
			if ((((search_it = search_ita) != database.end()) && ((search_it->second).size() > 0)) ||
				(((search_it = search_itb) != defaults.end()) && ((search_it->second).size() > 0)))
			{
				if ((search_it->second).size() == 1)
				{
					std::string s = (search_it->second)[0];

					//проверка на значение-переключатель
					if (((search_it = switchers.find(UpperCase(s))) != switchers.end()) && ((search_it->second).size() > 0))
						s = search_it->second[0];

					if ((UpperCase(s) == "FALSE") || (UpperCase(s) == "NO") || (s.c_str()[0] == '0'))
						res = false;
					else
						res = true;

					boolRes = true;
				}
				else
				{
					info('e') << "parameter " << name << " is list (only scalar is available)" << std::endl;
					exit(-1);
				}
			}
			else
				SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};//get(...)

			   		 	  
		/// \brief Считывание пары: ((строка, целое число), строка) вида ((str1, int), string) из базы данных
		/// 
		/// Метод get() для считывания пары: ((строка, целое число), строка) из базы данных
		///
		/// \param[in] name константная ссылка на строку --- имя считываемого параметра
		/// \param[out] res ссылка данные, считываемые из базы данных
		/// \param[in] defValue указатель на константу --- значение по умолчанию (по умолчанию nullptr)
		/// \param[in] echoDefault признак эхо-ответа при считывании значения по умолчанию (по умолчанию true)
		/// \return признак считывания переменной из базы данных (false - если считано значение по умолчанию)
		bool get(const std::string& name, std::pair<std::pair<std::string, int>, std::string>& res, const std::pair<std::pair<std::string, int>, std::string>* defValue = nullptr, bool echoDefault = true) const
		{
			bool boolRes = false;

			std::unordered_map<std::string, std::vector<std::string>>::const_iterator search_it, search_ita, search_itb;
			
			for (search_ita = database.begin(); search_ita != database.end(); ++search_ita)
				if (is_match(search_ita->first, UpperCase(name)))
					break;	

			for (search_itb = defaults.begin(); search_itb != defaults.end(); ++search_itb)
				if (is_match(search_itb->first, UpperCase(name)))
					break;
			
			//поиск ключа в базе и, если нет, то в базе умолчаний
			//if ((((search_it = database.find(UpperCase(name))) != database.end()) && ((search_it->second).size() > 0)) ||
			//	(((search_it = defaults.find(UpperCase(name))) != defaults.end()) && ((search_it->second).size() > 0)))			
			if ((((search_it = search_ita) != database.end()) && ((search_it->second).size() > 0)) ||
				(((search_it = search_itb) != defaults.end()) && ((search_it->second).size() > 0)))
			{
				if ((search_it->second).size() == 1)
				{
					std::string s = (search_it->second)[0];

					size_t posBegin = s.find('(');
					size_t posEnd = s.find(')');

					if ((posBegin != -1) && (posEnd == -1))
					{
						info('e') << "parameter " << name << " is given incorrectly" << std::endl;
						exit(-1);
					}
					else
					{
						std::string str1a, str1b, str2;
						str1a = s.substr(0, posBegin);

						if (((search_it = switchers.find(UpperCase(str1a))) != switchers.end()) && ((search_it->second).size() > 0))
							str1b = search_it->second[0];
						else
							str1b = "-1";

						std::stringstream ssStr1b;
						ssStr1b << str1b;
						int res1b;
						ssStr1b >> res1b;
									
						str2 = s.substr(posBegin + 1, posEnd - posBegin - 1);
					
						if ((posBegin == -1) && (posEnd == -1))
							res = { {str1a, res1b}, "" };
						else
							res = { {str1a, res1b}, str2 };

						boolRes = true;
					}
				}
				else
				{
					info('e') << "parameter " << name << " is list (only scalar is available)" << std::endl;
					exit(-1);
				}
			}
			else
				SetDefault(name, res, defValue, echoDefault);

			return boolRes;
		};
	};

}//namespace VMlib

#endif
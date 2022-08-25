/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.11   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2022/08/07     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Preprocessor.h                                                   |
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
\brief Заголовочный файл с описанием класса Preprocessor
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
*/ 

#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <string>
#include <sstream>

namespace VMlib
{

	/*!
	\brief Класс, позволяющий выполнять предварительную обработку файлов

	Выполняет следующие действия:

	- Исключает из файла однострочные и многостройчные комментарии в стиле С/С++
	- Разрывы строк заменяет на пробелы
	- Точки с запятой заменяет на разрывы строк

	\author Марчевский Илья Константинович
	\version 1.11
	\date 07 августа 2022 г.
	*/
	class Preprocessor
	{
	private:
		/// Признак встречи слэша внутри строки (внутри двойных кавычек)
		bool inStringLastEscape;

		/// Признак встречи слэша внутри символа (внутри одинарных кавычек)
		bool inCharLastEscape;

		/// Признак встречи слэша внутри однострочного комментария (после //)
		bool inInlineLastEscape;

		/// Признак встречи звездочки внутри многострочного комментария (после '/*')
		bool inMultilineLastStar;

		/// Опредедение currentParser как указателя на функцию-члена класса
		void(Preprocessor::*currentParser)(char ch, std::string& str);

		/// \brief Обработчик символа в режиме обычного состояния парcера
		///
		/// При появлении символов ", ', / переключает парсер в соответствующий режим
		/// Выводит обрабатываемый символ в поток (кроме слэша, который игнорируется)
		/// после его предварительной обработки фукцией Preprocessor::processSymbol
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void normalParser(char ch, std::string& str);

		/// \brief Обработчик символа в режиме парсера строки (внутри двойных кавычек)
		///
		/// При появлении парной кавычки переводит парсер в обычный режим, 
		/// обрабатываемый символ в неизменном виде передается в поток
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void inStringParser(char ch, std::string& str);

		/// \brief Обработчик символа в режиме парсера символа (внутри одинарных кавычек)
		///
		/// При появлении парного апострофа переводит парсер в обычный режим, 
		/// обрабатываемый символ в неизменном виде передается в поток
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void inCharParser(char ch, std::string& str);

		/// \brief Обработчик символа в режиме парсера выражения после слэша
		///
		/// При появлении символов / или * переключает парсер в режим обработки 
		/// однострочного или многострочного комментария
		/// \n Отальные символы выводит в поток после его предварительной 
		/// обработки фукцией Preprocessor::processSymbol
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void afterSlashParser(char ch, std::string& str);

		/// \brief Обработчик символа в режиме парсера однострочного комментария (после //)
		///
		/// Игнорирует все симводы до появления конца строки.
		/// \n При появлении конца строки переключает парсер в обычный режим
		/// и выводит конец строки в результирующую строчку
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void inInlineCommentParser(char ch, std::string& str);

		/// \brief Обработчик символа в режиме парсера многострочного комментария
		/// 
		/// Игнорирует все симводы до появления сочетания, соответствующего окончанию комментария.
		/// 
		/// При появлении символов окончания комментария выводит в результирующую строку пробел
		/// и переключает парсер в обычный режим
		/// \param[in] ch обрабатываемый символ
		/// \param[out] str ссылка на строку, в которую сохряняется результат
		void inMultilineCommentParser(char ch, std::string& str);

		/// \brief Базовая функция обработки символа
		/// В зависимости от входого символа возвращает на выход:
		/// - \n (конец строки) --- пробел
		/// - ; (точка с зарятой) --- конец строки
		/// - исходный символ в остальных случаях
		///
		/// \param[in] ch обрабатываемый символ
		/// return символ, который отправляется в обработанный файл
		char processSymbol(char ch);

	public:
		/// Строка, содержащая исходный файл в первоначальном виде
		std::string initialInput;

		/// Строка, содержащая результат промежуточной обработки файла
		std::string intermediateOutput;



		/// \brief Конструктор, принимающий на вход имя обрабатываемого файла
		///
		/// После посимвольной обработки формирует промежуточный вывод intermediateOutput
		/// \n Затем убирает из вывода пробелы и знаки табуляции и формирует 
		/// окончательный поток вывода resultStream
		/// \param[in] fileName константная ссылка на строку --- имя обрабатываемого файла
		Preprocessor(const std::string& fileName);

		/// Деструктор
		~Preprocessor() { };

		/// Строка, содержащая окончательный результат обработки файла
		std::string resultString;
		std::stringstream resultStream;

	};

}//namespace VMlib

#endif
/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.11   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2022/08/07     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: LogStream.h                                                      |
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
\brief Заголовочный файл с описанием класса LogStream
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
*/

#ifndef LOGSTREAM_H
#define LOGSTREAM_H

#include <ostream>
#include <string>

namespace VMlib
{

	/*!
	\brief Класс, определяющий работу с потоком логов
	\author Марчевский Илья Константинович
	\version 1.11
	\date 07 августа 2022 г.
	*/

	class LogStream
	{
	private:
		/// Указатель на поток вывода 
		std::ostream* pStr;

		/// Метка потока 
		std::string label;

	public:

		/// \brief Конструктор
		///
		/// Производит инициализацию потока логов пустым указателем
		LogStream()
			: pStr(nullptr), label("") {};


		/// \brief Связывание потока логов с потоком вывода
		///
		/// Устанавливает связь потока логов с каким-либо потоком вывода
		///
		/// \param[in] pStr_ указатель на поток вывода 
		/// \param[in] label_ константная ссылка на метку потока 
		void assignStream(std::ostream* pStr_, const std::string& label_)
		{
			pStr = pStr_;
			label = label_;
		}


		/// \brief Связывание потока логов с потоком вывода от другого потока логов
		///
		/// Устанавливает связь потока логов с потоком вывода, связанным с другим потоком логов
		///
		/// \param[in] infoStream_ существующий поток логов
		/// \param[in] label_ константная ссылка на дополнительную метку потока 
		void inheritStream(LogStream& infoStream_, const std::string& label_)
		{
			assignStream(infoStream_.pStr, infoStream_.label + "->" + label_);
		}

		/// Деструктор
		~LogStream() {};


		/// Вывод в поток логов пустой строки
		void endl()
		{
			*pStr << std::endl;
		}


		/// \brief Оператор вывода в поток логов
		///
		/// Оператор (), возврящающий ссылку на поток вывода и печатающий соответствующую метку в зависимости от признака
		/// Типы признаков:
		///   - i --- признак вывода информации, метка info
		///   - e --- признак вывода ошибки, метка ERROR
		///   - t --- признак вывода телеметрии, метка tele
		///   - - --- признак вывода без метки (два пробела в начале строки)
		///   - _ --- признак вывода без метки (четыре пробела в начале строки)
		///
		/// \param[in] c признак типа выводимой информации
		/// \return ссылку на поток вывода
		std::ostream& operator() (char c)
		{
			if (pStr != nullptr)
			{
				switch (c)
				{
				case 'i':
					*pStr << label << (label.length() > 0 ? " " : "") << "info: ";
					break;

				case 'e':
					*pStr << label << (label.length() > 0 ? " " : "") << "ERROR: ";
					break;

				case 't':
					*pStr << label << (label.length() > 0 ? " " : "") << "tele: ";
					break;

				case '-':
					*pStr << "  ";
					break;

				case '_':
					*pStr << "    ";
					break;

				}
				return *pStr;
			}
			else
			{
				pStr = new std::ostream(NULL);
				return *pStr;
			}
		}
	};

}//namespace VMlib;

#endif

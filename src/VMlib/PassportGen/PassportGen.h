/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.6    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/10/28     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: PassportGen.h                                                    |
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
\brief Заголовочный файл с описанием класса PasportGen
\author Марчевский Илья Константинович
\version 1.6
\date 28 октября 2019 г.
*/

#ifndef PASSPORTGEN_H
#define PASSPORTGEN_H

#include <iostream>
#include <memory>
#include <vector>

#include "LogStream.h"

namespace VMlib
{
	/*!
	\brief Структура, задающая параметры процесса интегрирования по времени
	\author Марчевский Илья Константинович
	\version 1.6
	\date 28 октября 2019 г.
	*/
	struct TimeDiscretizationProperties
	{
		/// Текущее время
		mutable double currTime;

		/// Начальное время 
		double timeStart;

		/// Конечное время
		double timeStop;

		/// Шаг по времени
		double dt;

		/// Число разрядов в имени файла
		int nameLength;

		/// Шаг сохранения кадров в текстовые файлы
		int saveTXT;

		/// Шаг сохранения кадров в бинарные файлы	
		int saveVTK;

		/// Шаг вычисления и сохранения скорости и давления
		int saveVP;
	};//TimeDiscretizationProperties


	/*!
	\brief Абстрактный класс, опеделяющий паспорт задачи
	\author Марчевский Илья Константинович
	\version 1.6
	\date 28 октября 2019 г.
	*/
	class PassportGen
	{
	protected:
		/// \brief Считывание всех параметров расчета из соответствующих потоков
		/// 
		/// \param[in] mainStream ссылка на основной поток
		/// \param[in] mechanicsStream ссылка на поток со словарем механических систем
		/// \param[in] defaultStream ссылка на поток с параметрами по умолчанию
		/// \param[in] switcherStream ссылка на поток со значениями параметров-переключателей
		/// \param[in] varsStream ссылка на поток с параметрами конкретной задачи и переменными
		virtual void GetAllParamsFromParser
		(
			std::istream& mainStream,
			std::istream& mechanicsStream,
			std::istream& defaultStream,
			std::istream& switcherStream,
			std::istream& varsStream
		) = 0;

		/// Поток для вывода логов и сообщений об ошибках
		mutable LogStream info;

		/// Печать всех параметров расчета в поток логов
		virtual void PrintAllParams() = 0;

	public:
		/// Рабочий каталог задачи
		std::string dir;

		/// Название задачи
		std::string problemName;
		
		/// Номер задачи
		size_t problemNumber;

		/// Структура с параметрами процесса интегрирования по времени
		TimeDiscretizationProperties timeDiscretizationProperties;

		/// \brief Конструктор
		///
		/// Осуществляет чтение всех данных из соответствующих потоков, полностью инициализирует паспорт
		///
		/// \param[in, out] infoStream базовый поток для вывода логов
		/// \param[in] _problemName константная ссылка наназвание задачи
		/// \param[in] _problemNumber номер (по счету) решаемой задачи
		/// \param[in] _filePassport константная ссылка на файл (без пути) с паспортом задачи
		/// \param[in] _mechanics константная ссылка на файл (c путем) со словарем механических систем
		/// \param[in] _defaults константная ссылка на имя файла (с путем) с параметрами по умолчанию
		/// \param[in] _switchers константная ссылка на имя файла (с путем) со значениями параметров-переключателей
		/// \param[in] vars константная ссылка на список переменных, заданных в виде строк
		PassportGen
		(
			LogStream& infoStream,
			const std::string& _problemName,
			const size_t _problemNumber,
			const std::string& _filePassport,
			const std::string& _mechanics,
			const std::string& _defaults,
			const std::string& _switchers,
			const std::vector<std::string>& vars
		);

		/// Деструктор
		virtual ~PassportGen() { };
	};

}//namespace VMlib


#endif


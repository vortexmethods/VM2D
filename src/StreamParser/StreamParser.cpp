/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: StreamParser.cpp                                                 |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.	                                                          |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса StreamParser
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "StreamParser.h"

//Создание экземпляров статических указателей
std::ostream* StreamParser::pinfo = defaults::defaultPinfo;
std::ostream* StreamParser::perr = defaults::defaultPerr;


//Конструктор, принимающий четыре потока
StreamParser::StreamParser(std::istream& mainStream, std::istream& defaultsStream, std::istream& switchersStream, std::istream& varsStream)
{
	ParseStream(mainStream, database);
	ParseStream(defaultsStream, defaults);
	ParseStream(switchersStream, switchers);
	ParseStream(varsStream, vars);
}//StreamParser(...)


//Конструктор, принимающий два потока
StreamParser::StreamParser(std::istream& mainStream, std::istream& defaultsStream)
{
	ParseStream(mainStream, database);
	ParseStream(defaultsStream, defaults);

	switchers.clear();
	vars.clear();
}//StreamParser(...)


//Конструктор, принимающий один поток
StreamParser::StreamParser(std::istream& mainStream, char openBracket, char closeBracket)
{
	ParseStream(mainStream, database, openBracket, closeBracket);
		
	defaults.clear();
	switchers.clear();
	vars.clear();
}//StreamParser(...)


//Pазбор строки, содержащей запятые, на отдельные строки
std::vector<std::string> StreamParser::StringToVector(std::string line, char openBracket, char closeBracket) //строка передается по значению, т.к. изменяется изнутри!
{
	if (line[0] != '{')
	{
		if (line.length() > 0)
			return { line };
		else
		{
			*pinfo << "parser info: WITH EMPTY PARAMETER" << std::endl;
			return std::vector<std::string>({});
		}
	}
	else
	{
		if (line[1] != '}')
		{
			std::vector<std::string> vecLine;
			size_t pos;
			while ( ((pos = line.find(',', 0)) < line.find(openBracket, 1)) ||
				    ((pos = line.find(',', std::max(1, (int)line.find(closeBracket, 1)))) != std::string::npos) ||
			      	(pos = line.length()-1)  )
			{
				std::string subline = line.substr(1, pos - 1);
				vecLine.push_back(subline);

				line.erase(1, subline.size() + 1);
				//info << "subline = !" << subline << "!" << endl;
			}
			return vecLine;
		}
		else
		{
			*pinfo << "parser info: WITH EMPTY LIST" << std::endl;
			return std::vector<std::string>({});
		}//else if (line[1] != '}')
	}//else if (st[0] != '{')
}//StringToVector(...)


//Объединение вектора (списка) из строк в одну строку
std::string StreamParser::VectorStringToString(const std::vector<std::string>& _vecString)
{
	std::string buf;
	for (size_t q = 0; q < _vecString.size(); ++q)
	{
		buf.append(_vecString[q]);
		buf.append("\n");
	}

	return buf;
}//VectorStringToString(...)


//Разбор строки на пару ключ-значение
std::pair<std::string, std::string> StreamParser::SplitString(std::string line)
{
	int posBegin = line.find('(', 0);
	int posEnd = line.find(')', 0);
	if ((posBegin != -1) && (posEnd == -1))
	{
		*perr << "parser error: ERROR WHILE PARSING LINE " << line << std::endl;
		exit(-1);
	}
	else
		return std::pair<std::string, std::string>(line.substr(0, posBegin), "{" + line.substr(posBegin + 1, posEnd - posBegin - 1) + "}");
}//SplitString(...)


//Парсинг заданного потока
void StreamParser::ParseStream(std::istream& streamName, std::unordered_map<std::string, std::vector<std::string>>& database, char openBracket, char closeBracket)
{
	std::string line, readline;

	while (streamName.good())
	{
		line = "";

		do
		{
			getline(streamName, readline);

			size_t posComment = std::min(readline.find('#'), readline.find("//"));
			//min, так как если не найдено -- то find возвращает string::npos, которое равно (-1), но приводится к большому положительному size_t
			if (posComment != std::string::npos)
				readline = readline.substr(0, posComment);
			
			line += readline;

			readline.erase(remove(readline.begin(), readline.end(), ' '), readline.end());
			readline.erase(remove(readline.begin(), readline.end(), '\t'), readline.end());

		} while ( (readline.size()>0) && (readline[readline.size() - 1] == '\\') && (streamName.good()));

		line.erase(remove(line.begin(), line.end(), '\\'), line.end());
		line.erase(remove(line.begin(), line.end(), ' '), line.end());
		line.erase(remove(line.begin(), line.end(), '\t'), line.end());

		int posEqual = line.find('=');
		if (posEqual != -1)
		{
			std::string key = line.substr(0, posEqual);
			auto search_it = database.find(key);

			if (search_it == database.end())
			{
				std::string value = line.substr(posEqual + 1, line.length());
				database.insert(make_pair(key, StringToVector(value, openBracket, closeBracket)));
			}
			else
			{
				*perr << "ERROR: KEY <" << key << "> FOUND TWICE" << std::endl;
				exit(-1);
			}
		}//if (posEqual != -1)
	}//while (streamName.good())

	//Указатель возвращаем в начало потока
	streamName.clear(); streamName.seekg(0, std::ios::beg);
}//ParseStream(...)

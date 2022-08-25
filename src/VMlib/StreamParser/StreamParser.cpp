/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.11   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2022/08/07     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: StreamParser.cpp                                                 |
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
\brief Файл кода с описанием класса StreamParser
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
*/

#include "StreamParser.h"

using namespace VMlib;

//Конструктор, принимающий четыре потока
StreamParser::StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, std::istream& defaultsStream, std::istream& switchersStream, std::istream& varsStream, std::vector<std::string> specificKey)
{
	info.inheritStream(infoStream, label);

	ParseStream(varsStream, vars);
	ParseStream(defaultsStream, defaults);
	ParseStream(switchersStream, switchers);
	ParseStream(mainStream, database, specificKey, true);	
}//StreamParser(...)


//Конструктор, принимающий два потока
StreamParser::StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, std::istream& defaultsStream, std::vector<std::string> specificKey)
{
	info.inheritStream(infoStream, label);

	ParseStream(mainStream, database, specificKey);
	ParseStream(defaultsStream, defaults);

	switchers.clear();
	vars.clear();
}//StreamParser(...)


//Конструктор, принимающий один поток
StreamParser::StreamParser(LogStream& infoStream, const std::string& label, std::istream& mainStream, char openBracket, char closeBracket)
{
	info.inheritStream(infoStream, label);

	ParseStream(mainStream, database, {}, false, openBracket, closeBracket);
		
	defaults.clear();
	switchers.clear();
	vars.clear();
}//StreamParser(...)

//Перевод строки в верхний регистр
std::string StreamParser::UpperCase(const std::string& line)
{
	std::string str(line);	
	bool flag = true;
	
	for (auto & c : str)
	{
		if (c == '"') 
			flag = !flag;
		if (flag)
			c = toupper(c);
	}

	//std::transform(str.begin(), str.end(), str.begin(), ::toupper);
	return str;
};//UpperCase(...)


//Pазбор строки, содержащей запятые, на отдельные строки
std::vector<std::string> StreamParser::StringToVector(std::string line, char openBracket, char closeBracket) //строка передается по значению, т.к. изменяется изнутри!
{
	//LogStream info;
	//info.inheritStream(infoStream, label);
	
	if (line[0] != '{')
	{
		if (line.length() > 0)
			return { line };
		else
		{
			//info('i') << "empty parameter" << std::endl;
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
				( (pos = line.find(',', std::max(1, static_cast<int>(line.find(closeBracket, 1)) ))) != std::string::npos) ||
			      	(pos = line.length()-1)  )
			{
				std::string subline = line.substr(1, pos - 1);
				if (subline != "")
					vecLine.push_back(subline);

				line.erase(1, subline.size() + 1);
			}
			return vecLine;
		}
		else
		{
			//info('i') << "empty parameter list" << std::endl;
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
std::pair<std::string, std::string> StreamParser::SplitString(LogStream& info, std::string line, bool upcase)
{
	size_t posBegin = line.find('(', 0);
	size_t posEnd = line.find(')', 0);
	if ((posBegin != -1) && (posEnd == -1))
	{
		info('e') << "error while parsing line " << line << std::endl;
		exit(-1);
	}
	else
	{
		std::string sub = line.substr(0, posBegin);
		if ((posBegin == -1) && (posEnd == -1))
			return std::pair<std::string, std::string>(upcase ? UpperCase(sub) : sub, "{}");
		else
			return std::pair<std::string, std::string>(upcase ? UpperCase(sub) : sub, "{" + line.substr(posBegin + 1, posEnd - posBegin - 1) + "}");
	}
}//SplitString(...)


//Парсинг заданного потока
void StreamParser::ParseStream(std::istream& stream, std::unordered_map<std::string, std::vector<std::string>>& database, std::vector<std::string> specificKey, bool replaceVars, char openBracket, char closeBracket)
{
	std::string line, readline;
	size_t defVarNum = 0;
		
	while (stream.good())
	{
		line = "";

		do
		{
			getline(stream, readline);

			size_t posComment = std::min(readline.find('#'), readline.find("//"));
			//min, так как если не найдено -- то find возвращает string::npos, которое равно (-1), но приводится к большому положительному size_t
			if (posComment != std::string::npos)
				readline = readline.substr(0, posComment);
			
			line += readline;

			readline.erase(remove(readline.begin(), readline.end(), ' '), readline.end());
			readline.erase(remove(readline.begin(), readline.end(), '\t'), readline.end());

		} while ( (readline.size()>0) && (readline[readline.size() - 1] == '\\') && (stream.good()));

		line.erase(remove(line.begin(), line.end(), '\\'), line.end());
		line.erase(remove(line.begin(), line.end(), ' '), line.end());
		line.erase(remove(line.begin(), line.end(), '\t'), line.end());

		size_t posEqual = line.find('=');
		if (posEqual != -1)
		{
			std::string key = line.substr(0, posEqual);
			auto search_it = database.find(UpperCase(key));

			if (search_it == database.end())
			{
				if ((specificKey.size() == 0) || (std::find(specificKey.begin(),specificKey.end(),key) != specificKey.end()) )
				{
					std::string value = line.substr(posEqual + 1, line.length());
					if (replaceVars)
						ReplaceVarsInString(value);
					database.insert(make_pair(UpperCase(key), StringToVector(value, openBracket, closeBracket)));
				}
			}
			else
			{
				info('e') << "key <" << key << "> is found twice" << std::endl;
				exit(-1);
			}
		}//if (posEqual != -1)
		else if (line.size() > 0)
		{
			std::string value = line;
			if (replaceVars)
				ReplaceVarsInString(value);

			std::stringstream ssDef;
			ssDef << "_defVar_" << defVarNum;

			database.insert(make_pair(UpperCase(ssDef.str()), StringToVector(value, openBracket, closeBracket)));

			++defVarNum;
		}
	}//while (streamName.good())

	//Указатель возвращаем в начало потока
	stream.clear(); 
	stream.seekg(0, std::ios::beg);
}//ParseStream(...)


// Замена переменных в строке их значениями
void StreamParser::ReplaceVarsInString(std::string& st)
{
	while (st.find("$") != -1)
	{
		size_t startVar = st.find('$');

		size_t endVar = -1 + std::min({
			st.find(" ", startVar + 1),
			st.find(",", startVar + 1),
			st.find(";", startVar + 1),
			st.find("$", startVar + 1),
			st.find("(", startVar + 1),
			st.find(")", startVar + 1),
			st.find("{", startVar + 1),
			st.find("}", startVar + 1),
			st.find("\n", startVar + 1),
			st.length()
		});

		std::string var = st.substr(startVar + 1, endVar - startVar);

		auto search_var = vars.find(UpperCase(var));
		if ((search_var != vars.end()) && ((search_var->second).size() > 0))
		{
			std::vector<std::string> findString = search_var->second;

			std::stringstream ss;

			ss << st.substr(0, startVar);
			
			if (findString.size() == 1)
			{
				ss << findString[0];
			}
			else
			{
				ss << "{" << findString[0];
				if (findString.size() > 0)
				{
					for (size_t sz = 1; sz < findString.size(); ++sz)
						ss << "," << findString[sz];
				}
				ss << "}";
			}


			size_t startP = endVar+1;
			size_t endP = st.length() - 1;

			if (startP <= endP)
				ss << st.substr(startP, endP);

			st = ss.str();
		}
		else
		{
			info('e') << "substitution $" << var << " is us undefined" << std::endl;
			exit(1);
		}
	}
}//ReplaceVarsInDataBaseValues(...)
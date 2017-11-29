/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: Preprocessor.cpp                                                 |
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
\brief Файл кода с описанием класса Preprocessor
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Preprocessor.h"

// Конструктор
Preprocessor::Preprocessor(const std::string& fileName)
	:
	inStringLastEscape(false),
	inCharLastEscape(false),
	inInlineLastEscape(false),
	inMultilineLastStar(false),
	currentParser(&Preprocessor::normalParser)
{
	FILE * inputFile = fopen(fileName.c_str(), "r");

	int symbol;

	while ((symbol = fgetc(inputFile)) != EOF)
	{
		initialInput.push_back(symbol);
		(this->*currentParser)(symbol, intermediateOutput);
	}

	fclose(inputFile);

	std::stringstream ss(intermediateOutput);

	std::string readline;

	while (ss.good())
	{
		getline(ss, readline);

		readline.erase(remove(readline.begin(), readline.end(), ' '), readline.end());
		readline.erase(remove(readline.begin(), readline.end(), '\t'), readline.end());

		resultStream << readline << std::endl;
	}

	resultString = resultStream.str();
}//Preprocessor(...)


// Базовая функция обработки символа
char Preprocessor::processSymbol(char ch)
{
	switch (ch)
	{
	case '\n':
		return ' ';

	case ';':
		return '\n';

	default:
		return ch;
	}
}//processSymbol(...)


// Обработчик символа в режиме обычного состояния парcера
void Preprocessor::normalParser(char ch, std::string& str)
{
	switch (ch)
	{
	case '\"':
		currentParser = &Preprocessor::inStringParser;
		break;

	case '\'':
		currentParser = &Preprocessor::inCharParser;
		break;

	case '/':
		currentParser = &Preprocessor::afterSlashParser;
		return;
	}
	str.push_back(processSymbol(ch));
}//normalParser(...)


// Обработчик символа в режиме парсера строки (внутри двойных кавычек)
void Preprocessor::inStringParser(char ch, std::string& str)
{
	str.push_back(ch);
	if (ch == '\"' && !inStringLastEscape)
		currentParser = &Preprocessor::normalParser;
	inStringLastEscape = (ch == '\\') && !inStringLastEscape;
}//inStringParser(...)


// Обработчик символа в режиме парсера символа (внутри одинарных кавычек)
void Preprocessor::inCharParser(char ch, std::string& str)
{
	str.push_back(ch);
	if (ch == '\'' && !inCharLastEscape)
		currentParser = &Preprocessor::normalParser;
	inCharLastEscape = (ch == '\\') && !inCharLastEscape;
}//inCharParser(...)


// Обработчик символа в режиме парсера выражения после слэша
void Preprocessor::afterSlashParser(char ch, std::string& str)
{
	switch (ch)
	{
	case '/':
		currentParser = &Preprocessor::inInlineCommentParser;
		return;

	case '*':
		currentParser = &Preprocessor::inMultilineCommentParser;
		return;

	default:
		str.push_back(processSymbol(ch));
	}
}//afterSlashParser(...)


// Обработчик символа в режиме парсера однострочного комментария (после //)
void Preprocessor::inInlineCommentParser(char ch, std::string& str)
{
	if (ch == '\n' && !inInlineLastEscape)
	{
		str.push_back(processSymbol(ch));
		currentParser = &Preprocessor::normalParser;
	}
	inInlineLastEscape = (ch == '\\') && !inInlineLastEscape;
}//inInlineCommentParser(...)


// Обработчик символа в режиме парсера многотрочного комментария
void Preprocessor::inMultilineCommentParser(char ch, std::string& str)
{
	if (ch == '/' && inMultilineLastStar)
	{
		str.push_back(' ');
		currentParser = &Preprocessor::normalParser;
	}
	inMultilineLastStar = (ch == '*');
}//inMultilineCommentParser(...)
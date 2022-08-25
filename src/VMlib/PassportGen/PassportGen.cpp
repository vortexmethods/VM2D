/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.11   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2022/08/07     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
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
\brief Файл кода с описанием класса PassportGen
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
*/

#include "defs.h"
#include "PassportGen.h"
#include "Preprocessor.h"
#include "StreamParser.h"

using namespace VMlib;

//Конструктор
PassportGen::PassportGen(LogStream& infoStream, const std::string& _problemName, const size_t _problemNumber, const std::string& _filePassport, const std::string& _mechanics, const std::string& _defaults, const std::string& _switchers, const std::vector<std::string>& vars)
: problemName(_problemName), problemNumber(_problemNumber)
{		
	std::stringstream ss;
	ss << "#" << problemNumber << " (" << problemName <<") passport";
	info.inheritStream(infoStream, ss.str());
	
	dir = "./" + problemName + "/";
}//PassportGen(...)
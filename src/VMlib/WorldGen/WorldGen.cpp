/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.7    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/11/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: WorldGen.cpp                                                     |
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
\brief Файл кода с описанием класса WorldGen
\author Марчевский Илья Константинович
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "WorldGen.h"

#include "PassportGen.h"

using namespace VMlib;

//Конструктор
WorldGen::WorldGen(const VMlib::PassportGen& passport_, const VMlib::Parallel& parallel_) :
	passportGen(passport_), parallel(parallel_)
{ };

//Функция, возвращающая признак завершения счета в решаемой задаче
bool WorldGen::isFinished() const  //Проверка условия, что счет завершен
{
	return (passportGen.timeDiscretizationProperties.currTime >= passportGen.timeDiscretizationProperties.timeStop);
}//isFinished()
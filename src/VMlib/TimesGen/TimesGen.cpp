/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: TimesGen.cpp                                                     |
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
\brief Файл кода с описанием класса TimersGen
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "TimesGen.h"

#include "PassportGen.h"

#include "WorldGen.h"

using namespace VMlib;

TimersGen::TimersGen(const WorldGen& W_, std::vector<std::string> labels)
    : W(W_)
{
    timerLabelList = labels;
    for (const auto& s : labels)
        timer.insert({ s, std::make_unique<vmTimer>(s) });
}//TimersGen(...)

void TimersGen::start(const std::string& timerLabel)
{
    try
    {
        timer.at(timerLabel)->start();
    }
    catch (...)
    {
        std::cout << "No timer \"" << timerLabel << "\" pre-determined!" << std::endl;
        exit(-16);
    }
}//start(...)

void TimersGen::stop(const std::string& timerLabel)
{
    try
    {
        timer.at(timerLabel)->stop();
    }
    catch (...)
    {
        std::cout << "No timer \"" << timerLabel << "\" pre-determined!" << std::endl;
        exit(-16);
    }
}//stop(...)

void TimersGen::resetAll()
{
    for (auto& tmr_ : timer)
        tmr_.second->reset();
}//resetAll(...)

void TimersGen::GenerateStatHeader()
{
    std::stringstream ss;
    ss << std::setw(6) << "Step" << "" \
        << std::setw(9) << "Time" << "" \
        << std::setw(9) << "N" << "";
    
    for (const auto& k : timerLabelList)
        ss << std::setw(9) << k << "";

    ss << std::setw(9) << "Other" << "";

    std::stringstream timeStatFileName;
    timeStatFileName << W.getPassportGen().dir << "timestat";

    std::ofstream timeStatFile(timeStatFileName.str());
    VMlib::PrintLogoToTextFile(timeStatFile, timeStatFileName.str(), "Time statistics (in milliseconds)");

    VMlib::PrintHeaderToTextFile(timeStatFile, ss.str());

    timeStatFile.close();
    timeStatFile.clear();
}//GenerateStatHeader()

void TimersGen::GenerateStatString(size_t stepNo, double curTime, size_t N)
{
    std::ofstream timestatFile(W.getPassportGen().dir + "timestat", std::ios::app);

    std::stringstream ss;
    ss  << std::setw(9) << stepNo << "" \
        << std::setw(9) << curTime << "" \
        //<< std::setw(9) << N << "";
        << std::setw(9) << W.nVtxBeforeMerging << "";
    double tOther = 0.0;

    ss.precision(2);
    ss << std::fixed;

    for (const auto& k : timerLabelList)
    {
        double duration = timer[k]->duration<vmTimer::ms>();
        ss << std::setw(9) << duration << "";
        tOther += (k == "Step" ? duration : -duration);
    }
    ss << std::setw(9) << tOther;

    timestatFile << std::endl << ss.str();
}//GenerateStatString()
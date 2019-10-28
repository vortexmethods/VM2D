/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.6    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/10/28     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: defs.cpp                                                         |
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
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\version 1.6   
\date 28 октября 2019 г.
*/

#include "defs.h"

using namespace VMlib;

std::string Months[12] = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

//Формирование строки с текущем временем и датой
std::string VMlib::CurrentDataTime()
{
	std::time_t t = std::time(nullptr);

#pragma warning(push)
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
	std::tm tm = *std::localtime(&t);
#pragma warning(pop)

	std::stringstream dateTimeStringStream;
	//dateTimeStringStream << "| This file was created automatically " << std::put_time(&tm, "%d %B %Y") << " at " << std::put_time(&tm, "%H:%M:%S");

	std::stringstream st_day, st_hour, st_min, st_sec;
	if (tm.tm_mday < 10)
		st_day << "0" << tm.tm_mday;
	else
		st_day << tm.tm_mday;

	if (tm.tm_hour < 10)
		st_hour << "0" << tm.tm_hour;
	else
		st_hour << tm.tm_hour;

	if (tm.tm_min < 10)
		st_min << "0" << tm.tm_min;
	else
		st_min << tm.tm_min;

	if (tm.tm_sec < 10)
		st_sec << "0" << tm.tm_sec;
	else
		st_sec << tm.tm_sec;

	dateTimeStringStream \
		<< st_day.str() << " " << Months[tm.tm_mon] << " " << 1900 + tm.tm_year \
		<< " at " \
		<< st_hour.str() << ":" << st_min.str() << ":" << st_sec.str();

	//std::string dateTimeString;
	//dateTimeStringStream >> dateTimeString;

	return dateTimeStringStream.str();
}

//Печать универсального логотипа
void VMlib::PrintUniversalLogoToStream(std::ostream& str)
{	
	str <<
			"-------------------------------------------------------------------------------" << '\n' << \
			"##\\    ##\\  ##\\      ##\\   ######\\   #######\\         ##\\  ######\\   #######\\  " << '\n' << \
			"## |   ## | ###\\    ### | ##  __##\\  ##  __##\\       ## | ## ___##\\  ##  __##\\ " << '\n' << \
			"## |   ## | ####\\  #### | \\__/  ## | ## |  ## |     ##  /  \\_ / ## | ## |  ## |" << '\n' << \
			"\\##\\  ##  | ##\\##\\## ## |  ######  | ## |  ## |    ##  /    ##### /  ## |  ## |" << '\n' << \
			" \\##\\##  /  ## \\###  ## | ##  ____/  ## |  ## |   ##  /     \\___##\\  ## |  ## |" << '\n' << \
			"  \\###  /   ## |\\#  /## | ## |       ## |  ## |  ##  /    ##\\   ## | ## |  ## |" << '\n' << \
			"   \\#  /    ## | \\_/ ## | ########\\  #######  | ##  /     \\######  | #######  |" << '\n' << \
			"    \\_/     \\__|     \\__| \\________| \\_______/  \\__/       \\______/  \\_______/ " << '\n' << \
			"-------------------------------------------------------------------------------" << '\n';
}

//Передача в поток вывода шапки программы VM2D/VM3D
void VMlib::PrintLogoToStream(std::ostream& str)
{
#ifdef CODE2D
		str <<
			"/*--------------------------------*- VM2D -*-----------------*---------------*\\" << '\n' << \
			"| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |" << '\n' << \
			"| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |" << '\n' << \
			"| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*" << '\n' << \
			"|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |" << '\n' << \
			"|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |" << '\n' << \
			"*-----------------------------------------------------------------------------*" << '\n';
#endif

#ifdef CODE3D
		str <<
			"/*--------------------------------*- VM3D -*-----------------*---------------*\\" << '\n' << \
			"| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |" << '\n' << \
			"| ##  ## ### ### ##  ## ##  ##  |  VM3D: Vortex Method       | 2019/05/30     |" << '\n' << \
			"| ##  ## ## # ##    ##  ##  ##  |  for 3D Flow Simulation    *----------------*" << '\n' << \
			"|  ####  ##   ## ##  ## ##  ##  |  Open Source Code                           |" << '\n' << \
			"|   ##   ##   ##  ####  #####   |  https://www.github.com/vortexmethods/VM3D  |" << '\n' << \
			"|                                                                             |" << '\n' << \
			"| Copyright (C) 2019 Ilia Marchevsky, Georgy Shcheglov, Sergey Dergachev      |" << '\n' << \
			"*-----------------------------------------------------------------------------*" << '\n';
#endif
}


//Формирование заголовка файла программы VM2D/VM3D
void VMlib::PrintLogoToTextFile(std::ofstream& str, const std::string& fileName, const std::string& descr)
{
	PrintLogoToStream(str);

	str <<
		"| File name: " << fileName;
	for (size_t q = 0; q < 65 - fileName.length(); ++q)
		str << " ";
	str <<
		"|"  << '\n';

	str <<
		"| Info: " << descr;
	for (size_t q = 0; q < 70 - descr.length(); ++q)
		str << " ";
	str <<
		"|" << '\n';
	
	
	std::string dateTimeString = "| This file was created automatically " + CurrentDataTime();
	//dateTimeStringStream >> dateTimeString;

	str << dateTimeString;	
	for (size_t q = 0; q < 78 - dateTimeString.length(); ++q)
		str << " ";
	str <<
		"|" << '\n';

	str <<
		"\\*---------------------------------------------------------------------------*/" << '\n';
}//PrintLogoToTextFile(...)




//Формирование подзаголовка в текстовом файле вывода программы VM2D/VM3D
void VMlib::PrintHeaderToTextFile(std::ofstream& str, const std::string& header)
{
	str << '\n';
	str << "// " << header << '\n';
	str << "//";
	for (size_t q = 0; q < header.length()+1; ++q)
		str << '-';
}

//Сохранение матрицы в поток
void VMlib::SaveToStream(const Eigen::MatrixXd& matr, std::ostream& str)
{
	str.precision(16);
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение комплекснозначной матрицы в поток
void VMlib::SaveToStream(const Eigen::MatrixXcd& matr, std::ostream& str)
{
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение вектора в поток
void VMlib::SaveToStream(const Eigen::VectorXd& vec, std::ostream& str)
{
	str.precision(16);
	for (int i = 0; i < vec.size(); ++i)
		str << vec(i) << " ";
}//SaveToStream(...)


//Сохранение списка из двумерных векторов (точек) в поток
void VMlib::SaveToStream(const std::vector<Point2D>& vec, std::ostream& str)
{
	for (size_t i = 0; i < vec.size(); ++i)
		str << "{ " << vec[i][0] << " " << vec[i][1] << " } ";
}//SaveToStream(...)


//Ядро сглаживания (Монагана)
double VMlib::M4(double t)
{
	double t2 = t*t;
	double mt = t > 0 ? t : -t;

	return (mt > 2.0) ? 0.0 : \
		(mt > 1.0) ? 0.5*sqr(2.0 - mt)*(1.0 - mt) : 1.0 - 2.5*t2 + 1.5*t2*mt;
}//M4(...)
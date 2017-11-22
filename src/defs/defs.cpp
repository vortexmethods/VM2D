/*!
\file
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include <iomanip>

#include "defs.h"

//Формирование заголовка файла программы VM2D
void PrintLogoToTextFile(std::ofstream& str, const std::string& fileName, const std::string& descr)
{
	str <<
		"/*--------------------------------*- C++ -*------------------*---------------*\\" << '\n' << \
		"| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |" << '\n' << \
		"| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |" << '\n' << \
		"| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*" << '\n' << \
		"|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |" << '\n' << \
		"|   ##   ##   ## ###### #####   |  http://www.github.com/vortexmethods/VM2D   |" << '\n' << \
		"|                                                                             |" << '\n' << \
		"| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |" << '\n' << \
		"*-----------------------------------------------------------------------------*" << '\n';

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

	std::time_t t = std::time(nullptr);

#pragma warning(push)
#pragma warning(disable : 4996) //_CRT_SECURE_NO_WARNINGS
	std::tm tm = *std::localtime(&t);
#pragma warning(pop)

	std::stringstream dateTimeStringStream;
	dateTimeStringStream << "| This file was created automatically " << std::put_time(&tm, "%d %B %Y") << " at " << std::put_time(&tm, "%H:%M:%S");
	//std::string dateTimeString;
	//dateTimeStringStream >> dateTimeString;

	str << dateTimeStringStream.str();		
	for (size_t q = 0; q < 78 - dateTimeStringStream.str().length(); ++q)
		str << " ";
	str <<
		"|" << '\n';

	str <<
		"\\*---------------------------------------------------------------------------*/" << '\n';
}//PrintLogoToTextFile(...)


//Формирование подзаголовка в текстовом файле вывода программы VM2D
void PrintHeaderToTextFile(std::ofstream& str, const std::string& header)
{
	str << '\n';
	str << "// " << header << '\n';
	str << "//";
	for (size_t q = 0; q < header.length()+1; ++q)
		str << '-';
}

//Сохранение матрицы в поток
void SaveToStream(const Eigen::MatrixXd& matr, std::ostream& str)
{
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение комплекснозначной матрицы в поток
void SaveToStream(const Eigen::MatrixXcd& matr, std::ostream& str)
{
	for (int i = 0; i < matr.rows(); ++i)
	{
		for (int j = 0; j < matr.cols(); ++j)
			str << matr(i, j) << " ";
		str << std::endl;
	}
}//SaveToStream(...)


//Сохранение вектора в поток
void SaveToStream(const Eigen::VectorXd& vec, std::ostream& str)
{
	for (int i = 0; i < vec.size(); ++i)
		str << vec(i) << " ";
}//SaveToStream(...)


//Сохранение списка из двумерных векторов (точек) в поток
void SaveToStream(const std::vector<Point2D>& vec, std::ostream& str)
{
	for (size_t i = 0; i < vec.size(); ++i)
		str << "{ " << vec[i][0] << " " << vec[i][1] << " } ";
}//SaveToStream(...)


//Ядро сглаживания (Монагана)
double M4(double t)
{
	double t2 = t*t;
	double mt = t > 0 ? t : -t;

	return (mt > 2.0) ? 0.0 : \
		(mt > 1.0) ? 0.5*sqr(2.0 - mt)*(1.0 - mt) : 1.0 - 2.5*t2 + 1.5*t2*mt;
}//M4(...)
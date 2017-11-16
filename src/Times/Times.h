/*!
\file
\brief Заголовочный файл с описанием класса Times
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef TIMES_H
#define TIMES_H

#include "defs.h"

/*!
\brief Класс для сбора статистики времени исполнения основных шагов алгоритма и вывода ее в файл
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Times
{
private:
	const Passport& passport;

	/// Обнуление одного временного периода 
	/// \param[out] period промежуток времени, начало и конец которого будут обнулены
	void ToZero(timePeriod& period)
	{
		period.first = 0;
		period.second = 0;
	}//ToZero(...)

public:
	/// Начало и конец процесса выполнения шага целиком
	timePeriod timeWholeStep;

	/// Начало и конец процесса выделения памяти под матрицу и правую часть
	timePeriod timeReserveMemoryForMatrixAndRhs;

	/// Начало и конец процесса заполнения матрицы и формирования правой части
	timePeriod timeFillMatrixAndRhs;

	/// Начало и конец процесса решения системы линейных алгебраических уравнений
	timePeriod timeSolveLinearSystem;

	/// Начало и конец процесса вычисления скоростей вихрей
	timePeriod timeCalcVortexVelo;

	/// Начало и конец процесса вычисления нагрузок
	timePeriod timeGetHydroDynamForce;

	/// Начало и конец процесса перемещения вихрей
	timePeriod timeMoveVortexes;

	/// Начало и конец процесса контроля протыкания
	timePeriod timeCheckInside;

	/// Начало и конец процесса реструктуризации пелены
	timePeriod timeRestruct;

	/// Начало и конец процесса сохранения кадра в файл
	timePeriod timeSaveKadr;
	
	/// Конструктор
	Times(const Passport& passport_)
		: passport(passport_) {};

	/// Деструктор
	~Times() {};
	
	/// Генерация заголовка файла временной статистики
	void GenerateStatHeader()
	{
		std::stringstream timeStatFileName;
		timeStatFileName << passport.dir << "timestat" << ".txt";

		std::ofstream timeStatFile(timeStatFileName.str());
		PrintLogoToTextFile(timeStatFile, timeStatFileName.str(), "Time statistics");

		PrintHeaderToTextFile(timeStatFile, "step Time N tStep tMem tMatRhs tSolve tVelo tForce tMove tInside tRestr tSave");

		timeStatFile.close();
		timeStatFile.clear();
	}//GenerateStatHeader()



	/// Сохранение строки со статистикой в файл временной статистики
	/// \param[in] currentStep номер текущего шага
	/// \param[in] N число вихрей в пелене
	void GenerateStatString(int currentStep, int N)
	{
		std::ofstream timestatFile(passport.dir + "timestat.txt", std::ios::app);
		
		timestatFile << std::endl
			<< currentStep << "\t"
			<< passport.physicalProperties.getCurrTime() << "\t"
			<< N << "\t"
			<< dT(timeWholeStep) << "\t"
			<< dT(timeReserveMemoryForMatrixAndRhs) << "\t"
			<< dT(timeFillMatrixAndRhs) << "\t"
			<< dT(timeSolveLinearSystem) << "\t"			
			<< dT(timeCalcVortexVelo) << "\t"
			<< dT(timeGetHydroDynamForce) << "\t"
			<< dT(timeMoveVortexes) << "\t"
			<< dT(timeCheckInside) << "\t"					
			<< dT(timeRestruct) << "\t"
			<< dT(timeSaveKadr);



		timestatFile.close();
	}//GenerateStatString(...)


	/// Обнуление состояния временной статистики
	void ToZero()
	{
		ToZero(timeWholeStep);
		ToZero(timeReserveMemoryForMatrixAndRhs);
		ToZero(timeFillMatrixAndRhs);
		ToZero(timeSolveLinearSystem);
		ToZero(timeCalcVortexVelo);
		ToZero(timeGetHydroDynamForce);
		ToZero(timeMoveVortexes);
		ToZero(timeCheckInside);
		ToZero(timeRestruct);
		ToZero(timeSaveKadr);		
	}


	/// Вычисление разницы во времени для пары засечек в секундах
	/// \param[in] константная ссылка на пару засечек времени
	/// \return разницу в секундах
	static double dT(const timePeriod& t)
	{
		return (double)(t.second - t.first) / CLOCKS_PER_SEC;
	}//dT(...)
};

#endif
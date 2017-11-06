/*!
\file
\brief Файл кода с описанием класса Parallel
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Parallel.h"


// Распределение задач по процессорам
void Parallel::SplitMPI(int n) const
{
	len.clear();
	disp.clear();

	int nPerP = n / nProcWork;

	for (int s = 0; s < nProcWork - 1; ++s)
	{
		len.push_back(nPerP);
		disp.push_back(s*nPerP);
	}
	len.push_back(n - nPerP * (nProcWork - 1));
	disp.push_back(nPerP * (nProcWork - 1));
}//SplitMPI(...)
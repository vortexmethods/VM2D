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

//#include <iostream>


// Распределение задач по процессорам
void Parallel::SplitMPI(int n) const
{	
	totalLen = n;
	MPI_Bcast(&totalLen, 1, MPI_INT, 0, commWork);
	
	if (myidWork == 0)
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
	}
	
	MPI_Scatter(len.data(), 1, MPI_INT, &myLen, 1, MPI_INT, 0, commWork);
	MPI_Scatter(disp.data(), 1, MPI_INT, &myDisp, 1, MPI_INT, 0, commWork);
}//SplitMPI(...)


/// Рассылка всего массива распределения витков цикла по всем процессорам
void Parallel::BCastAllLenDisp() const
{
	if (myidWork != 0)
	{ 
		len.resize(nProcWork);
		disp.resize(nProcWork);
	}
	MPI_Bcast(len.data(),  nProcWork, MPI_INT, 0, commWork);
	MPI_Bcast(disp.data(), nProcWork, MPI_INT, 0, commWork);
}
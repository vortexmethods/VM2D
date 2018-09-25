/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Parallel.cpp                                                     |
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
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Parallel
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/


#include "Parallel.h"


// Распределение задач по процессорам
parProp Parallel::SplitMPI(size_t n, bool bcastAll) const
{	
	parProp par;
	
	par.totalLen = (int)n;
	MPI_Bcast(&par.totalLen, 1, MPI_INT, 0, commWork);
	
	if (myidWork == 0)
	{
		par.len.clear();
		par.disp.clear();

		int nPerP = (int)(n / nProcWork);

		for (int s = 0; s < nProcWork - 1; ++s)
		{
			par.len.push_back(nPerP);
			par.disp.push_back(s*nPerP);
		}
		
		par.len.push_back((int)n - nPerP * (nProcWork - 1));
		par.disp.push_back(nPerP * (nProcWork - 1));
	}

	MPI_Scatter(par.len.data(), 1, MPI_INT, &par.myLen, 1, MPI_INT, 0, commWork);
	MPI_Scatter(par.disp.data(), 1, MPI_INT, &par.myDisp, 1, MPI_INT, 0, commWork);
	
	if (bcastAll)
	{
	    if (myidWork != 0)
	    { 
		par.len.resize(nProcWork);
		par.disp.resize(nProcWork);
	    }
	    MPI_Bcast(par.len.data(),  nProcWork, MPI_INT, 0, commWork);
	    MPI_Bcast(par.disp.data(), nProcWork, MPI_INT, 0, commWork);
	}
	
	return par;
	
}//SplitMPI(...)



/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.7    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/11/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Parallel.cpp                                                     |
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
\brief Файл кода с описанием класса Parallel
\author Марчевский Илья Константинович
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "Parallel.h"

using namespace VMlib;

// Распределение задач по процессорам
parProp Parallel::SplitMPIone(size_t n, bool bcastAll) const
{
	parProp par;

	par.totalLen = static_cast<int>(n);
	MPI_Bcast(&par.totalLen, 1, MPI_INT, 0, commWork);

	if (myidWork == 0)
	{
		par.len.clear();
		par.disp.clear();

		par.len.push_back(static_cast<int>(n));
		par.disp.push_back(0);
		
		for (int s = 1; s < nProcWork; ++s)
		{
			par.len.push_back(0);
			par.disp.push_back(static_cast<int>(n-1));
		}
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
		MPI_Bcast(par.len.data(), nProcWork, MPI_INT, 0, commWork);
		MPI_Bcast(par.disp.data(), nProcWork, MPI_INT, 0, commWork);
	}

	return par;

}//SplitMPIone(...)


// Распределение задач по процессорам
parProp Parallel::SplitMPI(size_t n, bool bcastAll) const
{	
	parProp par;
	
	par.totalLen = static_cast<int>(n);
	MPI_Bcast(&par.totalLen, 1, MPI_INT, 0, commWork);
	
	if (myidWork == 0)
	{
		par.len.clear();
		par.disp.clear();

		int nPerP = static_cast<int>(n / nProcWork);

		for (int s = 0; s < nProcWork - 1; ++s)
		{
			par.len.push_back(nPerP);
			par.disp.push_back(s*nPerP);
		}
		
		par.len.push_back(static_cast<int>(n) - nPerP * (nProcWork - 1));
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



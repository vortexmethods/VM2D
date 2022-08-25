/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: WakeDataBase2D.cpp                                               |
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
\brief Файл кода с описанием класса WakeDataBase
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#include "WakeDataBase2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Считывание вихревого следа из файла 
void WakeDataBase::ReadFromFile(const std::string& dir, const std::string& fileName)
{
	std::string filename = dir + fileName;
	std::ifstream wakeFile;
	
	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream wakeFile(VMlib::Preprocessor(filename).resultString);
		
		VMlib::LogStream XXX;
		VMlib::StreamParser wakeParser(XXX, "vortex wake file parser", wakeFile);
				
		wakeParser.get("vtx", vtx);
	}
}//ReadFromFile(...)


void WakeDataBase::SaveKadrVtk(const std::string& filePrefix) const
{
	W.getTimestat().timeSaveKadr.first += omp_get_wtime();

	if ((W.getParallel().myidWork == 0) && (W.ifDivisible(W.getPassport().timeDiscretizationProperties.saveVTK)))
	{
		std::string fname = VMlib::fileNameStep(filePrefix, W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");

		std::ofstream outfile;
		size_t numberNonZero = 0;

		if (vtx.size() > 0)
			numberNonZero += vtx.size();
		else
			for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				numberNonZero += W.getAirfoil(q).getNumberOfPanels();
		VMlib::CreateDirectory(W.getPassport().dir, "snapshots");

		outfile.open(W.getPassport().dir + "snapshots/" + fname);

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "VM2D VTK result: " << (W.getPassport().dir + "snapshots/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		outfile << "POINTS " << numberNonZero << " float" << std::endl;


		if (vtx.size() > 0)
			for (auto& v : vtx)
			{
				const Point2D& r = v.r();




				outfile << r[0] << " " << r[1] << " " << "0.0" << std::endl;

			}//for v		
		else
			for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
				{
					const Point2D& r = W.getAirfoil(q).getR(s);
					outfile << r[0] << " " << r[1] << " " << "0.0" << std::endl;
				}


		outfile << "CELLS " << numberNonZero << " " << 2 * numberNonZero << std::endl;
		for (size_t i = 0; i < numberNonZero; ++i)
			outfile << "1 " << i << std::endl;

		outfile << "CELL_TYPES " << numberNonZero << std::endl;
		for (size_t i = 0; i < numberNonZero; ++i)
			outfile << "1" << std::endl;

		outfile << std::endl;
		outfile << "POINT_DATA " << numberNonZero << std::endl;
		outfile << "SCALARS Gamma float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;

		if (vtx.size() > 0)
			for (auto& v : vtx)
				outfile << v.g() << std::endl;
		else
			for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
				{
					//const Point2D& r = W.getAirfoil(q).getR(s);
					outfile << "0.0" << std::endl;
				}

		outfile.close();
	}

	W.getTimestat().timeSaveKadr.second += omp_get_wtime();
}//SaveKadrVtk()

//MPI-синхронизация вихревого следа
void WakeDataBase::WakeSynchronize()
{
	W.getTimestat().timeOther.first += omp_get_wtime();
	
	int nV;
	if (W.getParallel().myidWork == 0)
		nV = (int)vtx.size();
	MPI_Bcast(&nV, 1, MPI_INT, 0, W.getParallel().commWork);

	if (W.getParallel().myidWork > 0)
		vtx.resize(nV);

	MPI_Bcast(vtx.data(), nV, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	W.getTimestat().timeOther.second += omp_get_wtime();
}//WakeSinchronize()
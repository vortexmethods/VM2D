/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "WakeDataBase2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"
#include "Gmres2D.h"

using namespace VM2D;

//Считывание вихревого следа из файла 
void WakeDataBase::ReadFromFile(const std::string& dir, const std::string& fileName)
{
	std::string filename = dir + fileName;
	std::ifstream wakeFile, testFile;
	
	char firstChar = '/', secondChar = '*';

	if (fileExistTest(filename, W.getInfo(), true, { "txt", "TXT" }))
	{		
		testFile.open(filename);
		testFile >> firstChar;
		testFile >> secondChar;

		//std::cout << "CHARS: " << firstChar << " " << secondChar << std::endl;
		testFile.close();
	}


	//Считывание из словаря
	if (firstChar == '/' && secondChar == '*')
	{
		if (fileExistTest(filename, W.getInfo(), true, { "txt", "TXT" }))
		{
			std::stringstream wakeFile(VMlib::Preprocessor(filename).resultString);

			VMlib::LogStream XXX;
			VMlib::StreamParser wakeParser(XXX, "vortex wake file parser", wakeFile);

			wakeParser.get("vtx", vtx);

			for (auto& v : vtx)
				v.sigma() = W.getPassport().wakeDiscretizationProperties.sigma0;

		}
	}
	else
	{
		//Считывание из обычного текстового файла
		if (fileExistTest(filename, W.getInfo(), true, { "txt", "TXT" }))
		{
			wakeFile.open(filename);
			int nnn;
			wakeFile >> nnn;
			vtx.reserve(nnn);
			for (int i = 0; i < nnn; ++i)
			{
				Vortex2D v;
				wakeFile >> v.r()[0] >> v.r()[1] >> v.g();				
				v.sigma() = W.getPassport().wakeDiscretizationProperties.sigma0;
				vtx.push_back(v);
			}

			wakeFile.close();
		}
	}

}//ReadFromFile(...)


void WakeDataBase::SaveKadrVtk(const std::string& filePrefix) const
{
	W.getTimers().start("Save");
	
	if (W.ifDivisible(W.getPassport().timeDiscretizationProperties.saveVtxStep))
	{		
		std::ofstream outfile;
		size_t numberNonZero = 0;

		if (vtx.size() > 0)
			numberNonZero += vtx.size();
		else
			for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				numberNonZero += W.getAirfoil(q).getNumberOfPanels();
		VMlib::CreateDirectory(W.getPassport().dir, "snapshots");

		if (W.getPassport().timeDiscretizationProperties.fileTypeVtx.second == 0) //text format vtk
		{
			std::string fname = VMlib::fileNameStep(filePrefix, W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");
			outfile.open(W.getPassport().dir + "snapshots/" + fname);

			outfile << "# vtk DataFile Version 2.0\n";
			outfile << "VM2D VTK result: " << (W.getPassport().dir + "snapshots/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << '\n';
			outfile << "ASCII\n";
			outfile << "DATASET UNSTRUCTURED_GRID\n";
			outfile << "POINTS " << numberNonZero << " float\n";

			
			if (vtx.size() > 0)
				for (auto& v : vtx)
				{
					const Point2D& r = v.r();
					outfile << r[0] << " " << r[1] << " " << "0.0\n";
				}//for v		
			else			
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
					{
						const Point2D& r = W.getAirfoil(q).getR(s);
						outfile << r[0] << " " << r[1] << " " << "0.0\n";
					}

			outfile << "CELLS " << numberNonZero << " " << 2 * numberNonZero << '\n';
			for (size_t i = 0; i < numberNonZero; ++i)
				outfile << "1 " << i << '\n';

			outfile << "CELL_TYPES " << numberNonZero << '\n';
			for (size_t i = 0; i < numberNonZero; ++i)
				outfile << "1\n";

			outfile << '\n';
			outfile << "POINT_DATA " << numberNonZero << '\n';
			outfile << "SCALARS Gamma float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			if (vtx.size() > 0)
				for (auto& v : vtx)
					outfile << v.g() << '\n';
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)											
						outfile << "0.0\n";					
						
			outfile << "SCALARS Sigma float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			if (vtx.size() > 0)
				for (auto& v : vtx)
					outfile << v.sigma() << '\n';
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
						outfile << "0.0\n";

#ifdef TURB
			outfile << "SCALARS nut float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			if (vtx.size() > 0)
				for (auto& v : vtx)
					outfile << sqr(v.sigma() / 4.48) / W.getPassport().timeDiscretizationProperties.dt - W.getPassport().physicalProperties.nu << '\n';
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
						outfile << "0.0\n";
#endif


			outfile.close();
		}//if fileType = text
		else if (W.getPassport().timeDiscretizationProperties.fileTypeVtx.second == 1) //binary format vtk
		{
			//Тест способа хранения чисел
			uint16_t x = 0x0001;
			bool littleEndian = (*((uint8_t*)&x));
			const char eolnBIN[] = "\n";

			std::string fname = VMlib::fileNameStep(filePrefix, W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");
			outfile.open(W.getPassport().dir + "snapshots/" + fname, std::ios::out | std::ios::binary);

			outfile << "# vtk DataFile Version 3.0" << "\r\n" << "VM2D VTK result: " << (W.getPassport().dir + "snapshots/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << eolnBIN;
			outfile << "BINARY" << eolnBIN;
			outfile << "DATASET UNSTRUCTURED_GRID" << eolnBIN << "POINTS " << numberNonZero << " " << "float" << eolnBIN;

			if (vtx.size() > 0)
			{
				Eigen::VectorXf rData = Eigen::VectorXf::Zero(vtx.size() * 3);
				for (size_t i = 0; i < vtx.size(); ++i)
				{
					rData(3 * i) = (float)(vtx[i].r())[0];
					rData(3 * i + 1) = (float)(vtx[i].r())[1];
				}//for i

				if (littleEndian)
					for (int i = 0; i < vtx.size() * 3; ++i)
						VMlib::SwapEnd(rData(i));
				outfile.write(reinterpret_cast<char*>(rData.data()), vtx.size() * 3 * sizeof(float));					
			}
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				{
					Eigen::VectorXf rData = Eigen::VectorXf::Zero(W.getAirfoil(q).getNumberOfPanels() * 3);
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
					{
						rData(3 * s) = (float)(W.getAirfoil(q).getR(s))[0];
						rData(3 * s + 1) = (float)(W.getAirfoil(q).getR(s))[1];
					}//for i

					if (littleEndian)
						for (int i = 0; i < W.getAirfoil(q).getNumberOfPanels() * 3; ++i)
							VMlib::SwapEnd(rData(i));
					outfile.write(reinterpret_cast<char*>(rData.data()), W.getAirfoil(q).getNumberOfPanels() * 3 * sizeof(float));
				}

			// CELLS
			std::vector<int> cells(2 * numberNonZero);
			for (size_t i = 0; i < numberNonZero; ++i)
			{
				cells[2 * i] = 1;
				cells[2 * i + 1] = (int)i;
			}

			std::vector<int> cellsTypes;
			cellsTypes.resize(numberNonZero, 1);

			if (littleEndian)
			{
				for (int i = 0; i < numberNonZero * 2; ++i)
					VMlib::SwapEnd(cells[i]);

				for (int i = 0; i < numberNonZero; ++i)
					VMlib::SwapEnd(cellsTypes[i]);
			}

			outfile << eolnBIN << "CELLS " << numberNonZero << " " << numberNonZero * 2 << eolnBIN;
			outfile.write(reinterpret_cast<char*>(cells.data()), numberNonZero * 2 * sizeof(int));
			outfile << eolnBIN << "CELL_TYPES " << numberNonZero << eolnBIN;
			outfile.write(reinterpret_cast<char*>(cellsTypes.data()), numberNonZero * sizeof(int));

			//gammas
			outfile << eolnBIN << "POINT_DATA " << numberNonZero << eolnBIN;
			outfile << eolnBIN << "SCALARS Gamma " << "float" << " 1" << eolnBIN;
			outfile << "LOOKUP_TABLE default" << eolnBIN;			

			if (vtx.size() > 0)
			{
				Eigen::VectorXf pData = Eigen::VectorXf::Zero(numberNonZero);
				for (int s = 0; s < vtx.size(); ++s)
					pData(s) = (float)vtx[s].g();
				
				if (littleEndian)
					for (int i = 0; i < vtx.size(); ++i)
						VMlib::SwapEnd(pData(i));
				outfile.write(reinterpret_cast<char*>(pData.data()), vtx.size() * sizeof(float));
			}
			else			
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				{
					Eigen::VectorXf pData = Eigen::VectorXf::Zero(W.getAirfoil(q).getNumberOfPanels());
					for (int s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
						pData(s) = 0;

					if (littleEndian)
						for (int i = 0; i < W.getAirfoil(q).getNumberOfPanels(); ++i)
							VMlib::SwapEnd(pData(i));
					outfile.write(reinterpret_cast<char*>(pData.data()), W.getAirfoil(q).getNumberOfPanels() * sizeof(float));					
				}
			
			outfile << eolnBIN;


			//sigmas
			outfile << eolnBIN << "SCALARS Sigma " << "float" << " 1" << eolnBIN;
			outfile << "LOOKUP_TABLE default" << eolnBIN;

			if (vtx.size() > 0)
			{
				Eigen::VectorXf pData = Eigen::VectorXf::Zero(numberNonZero);
				for (int s = 0; s < vtx.size(); ++s)
					pData(s) = (float)vtx[s].sigma();

				if (littleEndian)
					for (int i = 0; i < vtx.size(); ++i)
						VMlib::SwapEnd(pData(i));
				outfile.write(reinterpret_cast<char*>(pData.data()), vtx.size() * sizeof(float));
			}
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				{
					Eigen::VectorXf pData = Eigen::VectorXf::Zero(W.getAirfoil(q).getNumberOfPanels());
					for (int s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
						pData(s) = 0;

					if (littleEndian)
						for (int i = 0; i < W.getAirfoil(q).getNumberOfPanels(); ++i)
							VMlib::SwapEnd(pData(i));
					outfile.write(reinterpret_cast<char*>(pData.data()), W.getAirfoil(q).getNumberOfPanels() * sizeof(float));
				}

			outfile << eolnBIN;

#ifdef TURB
			//nut
			outfile << eolnBIN << "SCALARS nut " << "float" << " 1" << eolnBIN;
			outfile << "LOOKUP_TABLE default" << eolnBIN;

			if (vtx.size() > 0)
			{
				Eigen::VectorXf pData = Eigen::VectorXf::Zero(numberNonZero);
				for (int s = 0; s < vtx.size(); ++s)
					pData(s) = (float)(sqr(vtx[s].sigma() / 4.48) / W.getPassport().timeDiscretizationProperties.dt - W.getPassport().physicalProperties.nu);;

				if (littleEndian)
					for (int i = 0; i < vtx.size(); ++i)
						VMlib::SwapEnd(pData(i));
				outfile.write(reinterpret_cast<char*>(pData.data()), vtx.size() * sizeof(float));
			}
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
				{
					Eigen::VectorXf pData = Eigen::VectorXf::Zero(W.getAirfoil(q).getNumberOfPanels());
					for (int s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
						pData(s) = 0;

					if (littleEndian)
						for (int i = 0; i < W.getAirfoil(q).getNumberOfPanels(); ++i)
							VMlib::SwapEnd(pData(i));
					outfile.write(reinterpret_cast<char*>(pData.data()), W.getAirfoil(q).getNumberOfPanels() * sizeof(float));
				}

			outfile << eolnBIN;
#endif


			outfile.close();
		}//if binary

		if (W.getPassport().timeDiscretizationProperties.fileTypeVtx.second == 2) //csv
		{
			std::string fname = VMlib::fileNameStep(filePrefix, W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "csv");
			outfile.open(W.getPassport().dir + "snapshots/" + fname);

			outfile << "point,x,y,G,S";
#ifdef TURB
			outfile << ",nut";
#endif
			outfile << '\n';

			int counter = 0;

			if (vtx.size() > 0)
				for (auto& v : vtx)
				{
					outfile << counter++ << "," << v.r()[0] << "," << v.r()[1] << "," << v.g() << "," << v.sigma();
#ifdef TURB
					outfile << "," << sqr(v.sigma() / 4.48) / W.getPassport().timeDiscretizationProperties.dt - W.getPassport().physicalProperties.nu;
#endif
					outfile << '\n';
				}
			else
				for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					for (size_t s = 0; s < W.getAirfoil(q).getNumberOfPanels(); ++s)
					{
						const Point2D& r = W.getAirfoil(q).getR(s);
						outfile << counter++ << "," << r[0] << "," << r[1] << "," << "0.0" << "," << "0.0";
#ifdef TURB
						outfile << "," << 0.0;
#endif	
						outfile << '\n';
					}
			outfile.close();
		}//if fileType = text

	}

	W.getTimers().stop("Save");
}//SaveKadrVtk()

#ifdef TURB
void WakeDataBase::SaveScalarFields(const std::string& filePrefix,
	const std::vector<double>& VIP,
	const std::vector<double>& OIP,
	const std::vector<double>& nut,
	const std::vector<Point2D>& grNu,
	const std::vector<double>& RIP,
	const std::vector<double>& DELTAG
) const
{
	//W.getTimestat().timeSaveKadr.first += omp_get_wtime();

	if (W.ifDivisible(W.getPassport().timeDiscretizationProperties.saveVtxStep))
	{
		std::ofstream outfile;
		size_t numberNonZero = vtx.size();

		VMlib::CreateDirectory(W.getPassport().dir, "snapshots");

		//if (W.getPassport().timeDiscretizationProperties.fileTypeVtx.second == 0) //text format vtk
		{
			std::string fname = VMlib::fileNameStep(filePrefix, W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");
			outfile.open(W.getPassport().dir + "snapshots/" + fname);

			outfile << "# vtk DataFile Version 2.0\n";
			outfile << "VM2D VTK result: " << (W.getPassport().dir + "snapshots/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << '\n';
			outfile << "ASCII\n";
			outfile << "DATASET UNSTRUCTURED_GRID\n";
			outfile << "POINTS " << numberNonZero << " float\n";

			for (auto& v : vtx)
			{
				const Point2D& r = v.r();
				outfile << r[0] << " " << r[1] << " " << "0.0\n";
			}//for v		


			outfile << "CELLS " << numberNonZero << " " << 2 * numberNonZero << '\n';
			for (size_t i = 0; i < numberNonZero; ++i)
				outfile << "1 " << i << '\n';

			outfile << "CELL_TYPES " << numberNonZero << '\n';
			for (size_t i = 0; i < numberNonZero; ++i)
				outfile << "1\n";

			outfile << '\n';
			outfile << "POINT_DATA " << numberNonZero << '\n';
			outfile << "SCALARS Gamma float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : vtx)
				outfile << v.g() << '\n';


			outfile << "SCALARS Sigma float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : W.getWake().vtx)
				outfile << v.sigma() << '\n';

			outfile << "SCALARS VIP float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : VIP)
				outfile << v << '\n';

			outfile << "SCALARS OIP float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : OIP)
				outfile << v << '\n';


			outfile << "SCALARS nut float 1\n";
			outfile << "LOOKUP_TABLE default\n";


			for (auto& v : nut)
				outfile << v << '\n';

			outfile << "SCALARS RIP float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : RIP)
				outfile << v << '\n';


			outfile << "SCALARS DELTAG float 1\n";
			outfile << "LOOKUP_TABLE default\n";

			for (auto& v : DELTAG)
				outfile << v << '\n';


			outfile << "VECTORS grNu float\n";
			for (auto& v : grNu)
				outfile << v[0] << " " << v[1] << " 0.0\n";

			outfile.close();
		}//if fileType = text
	}

	//W.getTimestat().timeSaveKadr.second += omp_get_wtime();
}//SaveScalarFields(...)
#endif
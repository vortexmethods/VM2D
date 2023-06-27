/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: MeasureVP2D.cpp                                                  |
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
\brief Файл кода с описанием класса MeasureVP
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#if defined(_WIN32)
#include <direct.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

#include "MeasureVP2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "Boundary2DConstLayerAver.h"

#include "Mechanics2D.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Конструктор
MeasureVP::MeasureVP(const World2D& W_)
	: W(W_)
{
	wakeVP.reset(new WakeDataBase(W_));
};//MeasureVP(...)

//Чтение точек, в которых нужно посчитать давление и скорость
void MeasureVP::ReadPointsFromFile(const std::string& dir)
{
	std::string filename = dir + "pointsVP";
	std::ifstream VPFile;

	if (fileExistTest(filename, W.getInfo(), { "txt", "TXT" }))
	{
		std::stringstream VPFile(VMlib::Preprocessor(filename).resultString);

		VMlib::StreamParser VPParser(W.getInfo(), "velocity & pressure parser", VPFile);
		
		VPParser.get("points", initialPoints);
		VPParser.get("history", historyPoints);

		for (auto& pt : initialPoints)
			pt = pt.rotated(-W.getPassport().rotateAngleVpPoints * PI / 180.0);
		
		for (auto& pt : historyPoints)
			pt = pt.rotated(-W.getPassport().rotateAngleVpPoints * PI / 180.0);


		if (historyPoints.size() > 0)
		{
			VMlib::CreateDirectory(dir, "velPres");
			std::string VPFileNameList = W.getPassport().dir + "velPres/listPoints";
			std::ofstream VPFileList(VPFileNameList.c_str());

			VMlib::PrintLogoToTextFile(VPFileList, VPFileNameList.c_str(), "List of points, where velocity and pressure are measured in csv-files");

			VMlib::PrintHeaderToTextFile(VPFileList, "No.     FileName      X     Y");

			for (size_t q = 0; q < historyPoints.size(); ++q)
			{
				VPFileList << '\n' << q << '\t' << VMlib::fileNameStep("VP-atPoint-", 2, q, "csv") << '\t' << historyPoints[q][0] << '\t' << historyPoints[q][1];

				std::string VPFileNameCsv;
				VPFileNameCsv = W.getPassport().dir + "velPres/" + VMlib::fileNameStep("VP-atPoint-", 2, q, "csv");

				std::ofstream VPFileCsv(VPFileNameCsv.c_str());

				if (!W.getPassport().calcCoefficients)
					VPFileCsv << "t,Vx,Vy,p" << std::endl;
				else
					VPFileCsv << "t,CVx,CVy,Cp" << std::endl;
				
				VPFileCsv.close();
				VPFileCsv.clear();
			}

			VPFileList.close();
			VPFileList.clear();
		}
	}
}//ReadPointsFromFile(...)

// Инициализация векторов для вычисления скоростей и давлений
void MeasureVP::Initialization()
{
	W.getTimestat().timeVP.first += omp_get_wtime();

	if ((W.getPassport().timeDiscretizationProperties.saveVPstep > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVPstep)))
	{
		W.getInfo('i') << "Preparing VP points" << std::endl;


		wakeVP->vtx.clear();
		wakeVP->vtx.resize(0);

		Vortex2D addvtx;

		velocity.clear();
		velocity.resize(0);

		pressure.clear();
		pressure.resize(0);

		domainRadius.clear();
		domainRadius.resize(0);


		//отключаем проверку того, что точки расположены вне профиля
		/*
		//определяем точки, которые не находятся внутри профилей
		bool inside;
		for (size_t i = 0; i < initialPoints.size(); ++i)
		{
			inside = false;
			for (size_t j = 0; j < W.getNumberOfAirfoil(); ++j)
			{
				if (W.getAirfoil(j).isInsideGabarits(initialPoints[i]) && W.getAirfoil(j).IsPointInAirfoil(initialPoints[i]))
				{
					inside = true;
					break;
				}
			}
			if (!inside)
			{
				addvtx.r() = initialPoints[i];
				wakeVP.vtx.push_back(addvtx);
			}
		}
		*/

		//добавляем все точки, а не только те, которые вне профиля
		for (size_t i = 0; i < initialPoints.size(); ++i)
		{
			addvtx.r() = initialPoints[i];
			wakeVP->vtx.push_back(addvtx);
		}//for i

		for (size_t i = 0; i < historyPoints.size(); ++i)
		{
			addvtx.r() = historyPoints[i];
			wakeVP->vtx.push_back(addvtx);
		}//for i

	}
	W.getTimestat().timeVP.second += omp_get_wtime();

}//Initialization()

//Расчет поля давления
void MeasureVP::CalcPressure()
{
	int addWSize = (int)wakeVP->vtx.size();	
	pressure.resize(addWSize);	

	Point2D dri;
	Point2D Vi;
	Point2D vi;

	const Point2D& V0 = W.getPassport().physicalProperties.V0();
	const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;
	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	double P0 = /*1*/ 0.5 * V0.length2();

	Point2D cPan;

#pragma warning (push)
#pragma warning (disable: 4101)
	double alpha;
	double dst2eps;
#pragma warning (pop)

#pragma omp parallel for default(none) private(alpha, dri, Vi, vi, dst2eps, cPan) shared(P0, dt, eps2, V0, addWSize, std::cout, IDPI) 
	for (int i = 0; i < addWSize; ++i)
	{
		pressure[i] = P0;
		pressure[i] -= 0.5 * velocity[i].length2(); //2

		const Point2D& pt = wakeVP->vtx[i].r();

		for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
		{
			dri = pt - W.getWake().vtx[j].r();
			Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + V0;

			dst2eps = VMlib::boundDenom(dri.length2(), eps2);

			vi = IDPI * W.getWake().vtx[j].g() / dst2eps * dri.kcross();
			pressure[i] += vi & Vi; //3
		}

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			for (size_t j = 0; j < W.getBoundary(bou).virtualWake.vtx.size(); ++j)
			{
				dri = pt - W.getBoundary(bou).virtualWake.vtx[j].r();

				dst2eps = VMlib::boundDenom(dri.length2(), eps2);

				Vi = W.getBoundary(bou).virtualWake.vecHalfGamma[j];
				vi = IDPI * W.getBoundary(bou).virtualWake.vtx[j].g() / dst2eps * dri.kcross();

				pressure[i] += vi & Vi; //4
			}

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		{
			const Point2D& rcm = W.getAirfoil(bou).rcm;


			/*
			if ((i==0) && (bou == 2))
			{
				for (int q = 0; q < W.getAirfoil(bou).possibleWays.size(); ++q)
				{
					std::cout << "Possible way #" << q << ": " << std::endl;
					for (size_t s = 0; s < W.getAirfoil(bou).possibleWays[q].size(); ++s)
						std::cout << W.getAirfoil(bou).possibleWays[q][s] << " ";
					std::cout << std::endl;
				}
				
				std::ofstream of(W.getPassport().dir + "ways2.txt");
				for (size_t q = 0; q < W.getAirfoil(bou).getNumberOfPanels(); ++q)
				{
					cPan = 0.5 * (W.getAirfoil(bou).getR(q) + W.getAirfoil(bou).getR(q + 1));

					alpha = 0.0;

					int way = W.getAirfoil(bou).wayToVertex[q];
					if (way == 0)
						of << q << " " << rcm[0] << " " << rcm[1] << " " << cPan[0] << " " << cPan[1] << std::endl;
					else
					{
						of << q << " " << rcm[0] << " " << rcm[1];
						for (int q = 0; q < W.getAirfoil(bou).possibleWays[way - 1].size() + 1; ++q)
						{
							Point2D start = ((q == 0) ? rcm : W.getAirfoil(bou).possibleWays[way - 1][q - 1]);
							Point2D finish = ((q == W.getAirfoil(bou).possibleWays[way - 1].size()) ? cPan : W.getAirfoil(bou).possibleWays[way - 1][q]);
							of << " " << finish[0] << " " << finish[1];
						}
						of << std::endl;
					}
				}				
				of.close();
			}
			*/


			for (size_t j = 0; j < W.getAirfoil(bou).getNumberOfPanels(); ++j)
			{
				cPan = 0.5 * (W.getAirfoil(bou).getR(j) + W.getAirfoil(bou).getR(j + 1));
				
				alpha = 0.0;
				
				int way = W.getAirfoil(bou).wayToVertex[j];
				if (way < 0)				
					W.getInfo('i') << "Pressure computation is incorrect, way inside airfoil is undefined" << std::endl;				

				if (way == 0)
					alpha = atan2((cPan - pt) ^ (rcm - pt), (cPan - pt) & (rcm - pt));
				else
				{					
					for (int q = 0; q < W.getAirfoil(bou).possibleWays[way - 1].size() + 1; ++q)
					{
						Point2D start = ((q == 0) ? rcm : W.getAirfoil(bou).possibleWays[way - 1][q-1]);
						Point2D finish = ((q == W.getAirfoil(bou).possibleWays[way - 1].size()) ? cPan : W.getAirfoil(bou).possibleWays[way - 1][q]);
						alpha += atan2((finish - pt) ^ (start - pt), (finish - pt) & (start - pt));
					}
				}

				/// \todo Пока используем только средние значения свободного слоя на панелях
				pressure[i] += IDPI * alpha *	W.getBoundary(bou).sheets.freeVortexSheet(j, 0) * W.getAirfoil(bou).len[j] / dt;

				pressure[i] -= IDPI * alpha * W.getAirfoil(bou).gammaThrough[j] / dt;
			}
		}//for bou
	}//for i
		
	for (size_t i = 0; i < pressure.size(); ++i)
		pressure[i] *= W.getPassport().physicalProperties.rho;

}//CalcPressure()


//Сохранение в файл вычисленных скоростей и давлений
void MeasureVP::SaveVP()
{
	W.getTimestat().timeSaveKadr.first += omp_get_wtime();

	double scaleV = 1.0, scaleP = 1.0;
	if (W.getPassport().calcCoefficients)
	{
		const double& vRef = W.getPassport().physicalProperties.vRef;
		scaleV = 1.0 / vRef;
		scaleP = 1.0 / (0.5 * W.getPassport().physicalProperties.rho * sqr(vRef));
	}


	std::ofstream outfile;

	VMlib::CreateDirectory(W.getPassport().dir, "velPres");

	if (W.getPassport().timeDiscretizationProperties.fileTypeVP.second == 0) //text format VTK
	{
		std::string fname = VMlib::fileNameStep("VelPres", W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");
		outfile.open(W.getPassport().dir + "velPres/" + fname);

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "VM2D VTK result: " << (W.getPassport().dir + "velPres/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		outfile << "POINTS " << wakeVP->vtx.size() << " float" << std::endl;

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			double xi = (wakeVP->vtx[i].r())[0];
			double yi = (wakeVP->vtx[i].r())[1];
			outfile << xi << " " << yi << " " << "0.0" << std::endl;
		}//for i

		outfile << "CELLS " << wakeVP->vtx.size() << " " << 2 * wakeVP->vtx.size() << std::endl;
		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
			outfile << "1 " << i << std::endl;

		outfile << "CELL_TYPES " << wakeVP->vtx.size() << std::endl;
		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
			outfile << "1" << std::endl;

		outfile << std::endl;
		outfile << "POINT_DATA " << wakeVP->vtx.size() << std::endl;

		if (!W.getPassport().calcCoefficients)
			outfile << "VECTORS V float" << std::endl;
		else
			outfile << "VECTORS CV float" << std::endl;

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			outfile << velocity[i][0] * scaleV << " " << velocity[i][1] * scaleV << " 0.0" << std::endl;
		}//for i

		outfile << std::endl;

		if (!W.getPassport().calcCoefficients)
			outfile << "SCALARS P float 1" << std::endl;
		else
			outfile << "SCALARS CP float 1" << std::endl;

		outfile << "LOOKUP_TABLE default" << std::endl;

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			outfile << pressure[i] * scaleP << std::endl;
		}//for i

		outfile.close();
	}
	else if (W.getPassport().timeDiscretizationProperties.fileTypeVP.second == 1) //binary format VTK
	{
		//Тест способа хранения чисел
		uint16_t x = 0x0001;
		bool littleEndian = (*((uint8_t*)&x));
		const char eolnBIN[] = "\n";

		std::string fname = VMlib::fileNameStep("VelPres", W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "vtk");
		outfile.open(W.getPassport().dir + "velPres/" + fname, std::ios::out | std::ios::binary);

		outfile << "# vtk DataFile Version 3.0" << "\r\n" << "VM2D VTK result: " << (W.getPassport().dir + "velPres/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << eolnBIN;
		outfile << "BINARY" << eolnBIN;
		outfile << "DATASET UNSTRUCTURED_GRID" << eolnBIN << "POINTS " << wakeVP->vtx.size() << " " << "float" << eolnBIN;


		Eigen::VectorXf pData = Eigen::VectorXf::Zero(wakeVP->vtx.size());
		Eigen::VectorXf vData = Eigen::VectorXf::Zero(wakeVP->vtx.size() * 3);
		Eigen::VectorXf rData = Eigen::VectorXf::Zero(wakeVP->vtx.size() * 3);


		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			rData(3 * i + 0) = (float)(wakeVP->vtx[i].r())[0];
			rData(3 * i + 1) = (float)(wakeVP->vtx[i].r())[1];
		}//for i

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			vData(3 * i + 0) = (float)(velocity[i][0] * scaleV);
			vData(3 * i + 1) = (float)(velocity[i][1] * scaleV);
		}//for i

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			pData(i) = (float)(pressure[i] * scaleP);
		}//for i

		//POINTS
		if (littleEndian)
			for (int i = 0; i < wakeVP->vtx.size() * 3; ++i)
				VMlib::SwapEnd(rData(i));
		outfile.write(reinterpret_cast<char*>(rData.data()), wakeVP->vtx.size() * 3 * sizeof(float));

		// CELLS
		std::vector<int> cells(2 * wakeVP->vtx.size());
		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			cells[2 * i + 0] = 1;
			cells[2 * i + 1] = (int)i;
		}

		std::vector<int> cellsTypes;
		cellsTypes.resize(wakeVP->vtx.size(), 1);

		if (littleEndian)
		{
			for (int i = 0; i < wakeVP->vtx.size() * 2; ++i)
				VMlib::SwapEnd(cells[i]);

			for (int i = 0; i < wakeVP->vtx.size(); ++i)
				VMlib::SwapEnd(cellsTypes[i]);
		}

		outfile << eolnBIN << "CELLS " << wakeVP->vtx.size() << " " << wakeVP->vtx.size() * 2 << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cells.data()), wakeVP->vtx.size() * 2 * sizeof(int));
		outfile << eolnBIN << "CELL_TYPES " << wakeVP->vtx.size() << eolnBIN;
		outfile.write(reinterpret_cast<char*>(cellsTypes.data()), wakeVP->vtx.size() * sizeof(int));

		//VECTORS V
		if (littleEndian)
			for (int i = 0; i < wakeVP->vtx.size() * 3; ++i)
				VMlib::SwapEnd(vData(i));

		outfile << eolnBIN << "POINT_DATA " << wakeVP->vtx.size() << eolnBIN;

		if (!W.getPassport().calcCoefficients)
			outfile << "VECTORS V " << "float" << eolnBIN;
		else
			outfile << "VECTORS CV " << "float" << eolnBIN;

		outfile.write(reinterpret_cast<char*>(vData.data()), wakeVP->vtx.size() * 3 * sizeof(float));

		//SCALARS P	
		if (littleEndian)
			for (int i = 0; i < wakeVP->vtx.size(); ++i)
				VMlib::SwapEnd(pData(i));

		if (!W.getPassport().calcCoefficients)
			outfile << eolnBIN << "SCALARS P " << "float" << " 1" << eolnBIN;
		else
			outfile << eolnBIN << "SCALARS CP " << "float" << " 1" << eolnBIN;

		outfile << "LOOKUP_TABLE default" << eolnBIN;
		outfile.write(reinterpret_cast<char*>(pData.data()), wakeVP->vtx.size() * sizeof(float));

		outfile << eolnBIN;

		outfile.close();
	}
	else if (W.getPassport().timeDiscretizationProperties.fileTypeVP.second == 2) //csv
	{
		std::string fname = VMlib::fileNameStep("VelPres", W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "csv");
		outfile.open(W.getPassport().dir + "velPres/" + fname);

		if (!W.getPassport().calcCoefficients)
			outfile << "point,x,y,Vx,Vy,P" << std::endl;
		else
			outfile << "point,x,y,CVx,CVy,CP" << std::endl;

		for (size_t i = 0; i < wakeVP->vtx.size(); ++i)
		{
			double xi = (wakeVP->vtx[i].r())[0];
			double yi = (wakeVP->vtx[i].r())[1];
			outfile << i << "," << xi << "," << yi \
				<< "," << velocity[i][0] * scaleV << "," << velocity[i][1] * scaleV \
				<< "," << pressure[i] * scaleP << std::endl;
		}//for i
		outfile.close();
	}
	else if (W.getPassport().timeDiscretizationProperties.fileTypeVP.second == 3) //csvBundle
	{
		std::string fnameBunCsv = W.getPassport().dir + "velPres/" + "velPresBundle.csv";
		std::ofstream VPFileBunCsv(fnameBunCsv.c_str(), W.getCurrentStep() ? std::ios::app : std::ios::out);
		{
			if (W.getCurrentStep() == 0)
				VPFileBunCsv << wakeVP->vtx.size() << std::endl;


			for (size_t q = 0; q < wakeVP->vtx.size(); ++q)
			{
				VPFileBunCsv << W.getCurrentStep() << "," \
					<< W.getPassport().physicalProperties.getCurrTime() << "," \
					<< q << "," \
					<< wakeVP->vtx[q].r()[0] << "," << wakeVP->vtx[q].r()[1] << "," \
					<< velocity[q][0] * scaleV << "," << velocity[q][1] * scaleV << "," \
					<< pressure[q] * scaleP << std::endl;
			}
		}

	}


	//Вывод в csv-файлы точек "historyPoints"
	for (size_t q = 0; q < historyPoints.size(); ++q)
	{
		std::string VPFileNameCsv;
		VPFileNameCsv = W.getPassport().dir + "velPres/" + VMlib::fileNameStep("VP-atPoint-", 2, q, "csv");

		std::ofstream VPFileCsv(VPFileNameCsv.c_str(), W.getCurrentStep() ? std::ios::app : std::ios::out);

		if (W.getCurrentStep() == 0)
		{
			if (!W.getPassport().calcCoefficients)
				VPFileCsv << "point,time,Vx,Vy,P" << std::endl;
			else
				VPFileCsv << "point,time,CVx,CVy,CP" << std::endl;
		}

		VPFileCsv << q << "," << W.getPassport().physicalProperties.getCurrTime() << ","
			<< velocity[initialPoints.size() + q][0] * scaleV << ","
			<< velocity[initialPoints.size() + q][1] * scaleV << ","
			<< pressure[initialPoints.size() + q] * scaleP << std::endl;
		VPFileCsv.close();
	}

	W.getTimestat().timeSaveKadr.second += omp_get_wtime();
}//SaveVP()
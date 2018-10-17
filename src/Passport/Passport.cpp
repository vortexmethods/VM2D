/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Passport.cpp                                                     |
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
\brief Файл кода с описанием класса Passport
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#include "defs.h"
#include "Passport.h"
#include "Preprocessor.h"
#include "StreamParser.h"

//Конструктор
Passport::Passport(LogStream& infoStream, const std::string& _problemName, const size_t _problemNumber, const std::string& _filePassport, const std::string& _mechanics, const std::string& _defaults, const std::string& _switchers, const std::vector<std::string>& vars)
: problemName(_problemName), problemNumber(_problemNumber)
{		
	std::stringstream ss;
	ss << "#" << problemNumber << " (" << problemName <<") passport";
	info.inheritStream(infoStream, ss.str());
	
	dir = "./" + problemName + "/";

	std::string fileFullName = dir + _filePassport;
	std::string mechanicsFileFullName = _mechanics;
	std::string defaultsFileFullName = _defaults;
	std::string switchersFileFullName = _switchers;
	
	if (fileExistTest(fileFullName, info) 
		&& fileExistTest(defaultsFileFullName, info)
		&& fileExistTest(switchersFileFullName, info)
		&& fileExistTest(mechanicsFileFullName, info))
	{
		std::string str = StreamParser::VectorStringToString(vars);
		std::stringstream varsStream(str);
		
		std::stringstream mainStream;
		mainStream << Preprocessor(fileFullName).resultString;
		
		std::stringstream mechanicsStream;
		mechanicsStream << Preprocessor(mechanicsFileFullName).resultString;

		std::stringstream defaultsStream;
 		defaultsStream << Preprocessor(defaultsFileFullName).resultString;
		
		std::stringstream switchersStream;
		switchersStream << Preprocessor(switchersFileFullName).resultString;

		GetAllParamsFromParser(mainStream, mechanicsStream, defaultsStream, switchersStream, varsStream);

		PrintAllParams();
	}	
}//Passport(...)


//Считывание всех параметров расчета из соответствующих потоков
void Passport::GetAllParamsFromParser
(
	std::istream& mainStream, 
	std::istream& mechanicsStream,
	std::istream& defaultStream, 
	std::istream& switcherStream, 
	std::istream& varsStream
)
{

// 1. Разбор паспорта в целом
	
	//создаем парсер и связываем его с нужным потоком
	std::unique_ptr<StreamParser> parser;
	parser.reset(new StreamParser(info, "parser", mainStream, defaultStream, switcherStream, varsStream));
		
	//считываем общие парамеры
	parser->get("rho", physicalProperties.rho);
	parser->get("vInf", physicalProperties.vInf);
	parser->get("timeAccel", physicalProperties.timeAccel, &defaults::defaultTimeAccel);
	parser->get("nu", physicalProperties.nu);

	parser->get("timeStart", timeDiscretizationProperties.timeStart, &defaults::defaultTimeStart);
	parser->get("timeStop", timeDiscretizationProperties.timeStop);
	parser->get("dt", timeDiscretizationProperties.dt);
	parser->get("saveTXT", timeDiscretizationProperties.saveTXT, &defaults::defaultSaveTXT);
	parser->get("saveVTK", timeDiscretizationProperties.saveVTK, &defaults::defaultSaveVTK);
	parser->get("saveVP",  timeDiscretizationProperties.saveVP,  &defaults::defaultSaveVP);


	parser->get("eps", wakeDiscretizationProperties.eps);
	wakeDiscretizationProperties.eps2 = wakeDiscretizationProperties.eps* wakeDiscretizationProperties.eps;
	parser->get("epscol", wakeDiscretizationProperties.epscol);
	parser->get("distFar", wakeDiscretizationProperties.distFar, &defaults::defaultDistFar);
	parser->get("delta", wakeDiscretizationProperties.delta, &defaults::defaultDelta);
	parser->get("vortexPerPanel", wakeDiscretizationProperties.vortexPerPanel, &defaults::defaultVortexPerPanel);
	
	parser->get("linearSystemSolver", numericalSchemes.linearSystemSolver);
	parser->get("velocityComputation", numericalSchemes.velocityComputation);
	parser->get("wakeMotionIntegrator", numericalSchemes.wakeMotionIntegrator);
	
	parser->get("airfoilsDir", airfoilsDir, &defaults::defaultAirfoilsDir);
	parser->get("wakesDir",    wakesDir,    &defaults::defaultWakesDir);

	parser->get("fileWake", wakeDiscretizationProperties.fileWake, &defaults::defaultFileWake);

	std::vector<std::string> airfoil;
	parser->get("airfoil", airfoil, &defaults::defaultAirfoil);

	// 2. Разбор параметров профилей

	//определяем число профилей и организуем цикл по ним
	size_t nAirfoil = airfoil.size();
	//*(defaults::defaultPinfo) << "Number of airfoils = " << nAirfoil << endl;
	for (size_t i = 0; i < nAirfoil; ++i) 
	{
		//делим имя файла + выражение в скобках на 2 подстроки
		std::pair<std::string, std::string> airfoilLine = StreamParser::SplitString(info, airfoil[i], false);

		AirfoilParams prm;
		//первая подстрока - имя файла
		prm.fileAirfoil = airfoilLine.first;

		//вторую подстроку разделяем на вектор из строк по запятым, стоящим вне фигурных скобок		
		std::vector<std::string> vecAirfoilLineSecond = StreamParser::StringToVector(airfoilLine.second, '{', '}');
				
		//создаем парсер и связываем его с параметрами профиля
		std::stringstream aflStream(StreamParser::VectorStringToString(vecAirfoilLineSecond));
		std::unique_ptr<StreamParser> parserAirfoil;

		parserAirfoil.reset(new StreamParser(info, "airfoil parser", aflStream, defaultStream, switcherStream, varsStream));

		//считываем нужные параметры с учетом default-значений
		parserAirfoil->get("basePoint", prm.basePoint, &defaults::defaultBasePoint);
		parserAirfoil->get("scale", prm.scale, &defaults::defaultScale);
		parserAirfoil->get("angle", prm.angle, &defaults::defaultAngle);
		prm.angle *= PI / 180.0;
		parserAirfoil->get("panelsType", prm.panelsType, &defaults::defaultPanelsType);
		parserAirfoil->get("boundaryConditionSatisfaction", prm.boundaryCondition, &defaults::defaultBoundaryCondition);
		parserAirfoil->get("mechanicalSystem", prm.mechanicalSystem, &defaults::defaultMechanicalSystem);

		if (prm.mechanicalSystem == defaults::defaultMechanicalSystem)
		{
			prm.mechanicalSystemType = 0;
			prm.mechanicalSystemParameters = "";
		}
		else
		{
			std::unique_ptr<StreamParser> parserMechanicsList;
			std::unique_ptr<StreamParser> parserSwitchers;
			parserMechanicsList.reset(new StreamParser(info, "mechanical parser", mechanicsStream, defaultStream, switcherStream, varsStream, { prm.mechanicalSystem }));
			parserSwitchers.reset(new StreamParser(info, "switchers parser" ,switcherStream));

			std::string mechString;

			parserMechanicsList->get(prm.mechanicalSystem, mechString);
			
			//делим тип мех.системы + выражение в скобках (ее параметры) на 2 подстроки
			std::pair<std::string, std::string> mechanicsLine = StreamParser::SplitString(info, mechString);

			std::string mechTypeAlias = mechanicsLine.first;
			parserSwitchers->get(mechTypeAlias, prm.mechanicalSystemType);

			//вторую подстроку разделяем на вектор из строк по запятым, стоящим вне фигурных скобок		
			std::vector<std::string> vecMechLineSecond = StreamParser::StringToVector(mechanicsLine.second, '{', '}');
			prm.mechanicalSystemParameters = StreamParser::VectorStringToString(vecMechLineSecond);
		}

		//отправляем считанные параметры профиля в структуру данных паспорта 
		airfoilParams.push_back(prm);

	} //for i
}//GetAllParamsFromParser(...)


//Печать всех параметров расчета в поток логов
void Passport::PrintAllParams()
{
	const std::string str = "passport info: ";
	
	info('i') << "--- Passport info ---" << std::endl;
	info('-') << "rho = " << physicalProperties.rho << std::endl;
	info('-') << "vInf = " << physicalProperties.vInf << std::endl;
	info('-') << "nu = " << physicalProperties.nu << std::endl;
	info('-') << "timeStart = " << timeDiscretizationProperties.timeStart << std::endl;
	info('-') << "timeStop = " << timeDiscretizationProperties.timeStop << std::endl;
	info('-') << "dt = " << timeDiscretizationProperties.dt << std::endl;
	info('-') << "saveTXT = " << timeDiscretizationProperties.saveTXT << std::endl;
	info('-') << "saveVTK = " << timeDiscretizationProperties.saveVTK << std::endl;
	info('-') << "saveVP = " << timeDiscretizationProperties.saveVP << std::endl;
	info('-') << "eps = " << wakeDiscretizationProperties.eps << std::endl;
	info('-') << "eps2 = " << wakeDiscretizationProperties.eps2 << std::endl;
	info('-') << "epscol = " << wakeDiscretizationProperties.epscol << std::endl;
	info('-') << "distFar = " << wakeDiscretizationProperties.distFar << std::endl;
	info('-') << "delta = " << wakeDiscretizationProperties.delta << std::endl;
	info('-') << "vortexPerPanel = " << wakeDiscretizationProperties.vortexPerPanel << std::endl;
	info('-') << "linearSystemSolver = " << numericalSchemes.linearSystemSolver << std::endl;
	info('-') << "velocityComputation = " << numericalSchemes.velocityComputation << std::endl;
	info('-') << "wakeMotionIntegrator = " << numericalSchemes.wakeMotionIntegrator << std::endl;
	
	info('-') << "airfoilsDir = " << airfoilsDir << std::endl;
	info('-') << "wakesDir = " << wakesDir << std::endl;

	info('-') << "number of airfoils = " << airfoilParams.size() << std::endl;
	for (size_t q = 0; q < airfoilParams.size(); ++q)
	{
		info('_') << "airfoil[" << q << "]_file = " << airfoilParams[q].fileAirfoil << std::endl;
		info('_') << "airfoil[" << q << "]_basePoint = " << airfoilParams[q].basePoint << std::endl;
		info('_') << "airfoil[" << q << "]_scale = " << airfoilParams[q].scale << std::endl;
		info('_') << "airfoil[" << q << "]_angle = " << airfoilParams[q].angle << std::endl;
		info('_') << "airfoil[" << q << "]_panelType = " << airfoilParams[q].panelsType << std::endl;
		info('_') << "airfoil[" << q << "]_boundaryCondition = " << airfoilParams[q].boundaryCondition << std::endl;
		info('_') << "airfoil[" << q << "]_mechanicalSystem = " << airfoilParams[q].mechanicalSystem << std::endl;
	}

	info('-') << "fileWake = " << wakeDiscretizationProperties.fileWake << std::endl;
}//PrintAllParams()
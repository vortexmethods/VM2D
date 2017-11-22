/*!
\file
\brief Файл кода с описанием класса Passport
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

//#include <stringstream>

#include "Passport.h"

//Конструктор
Passport::Passport(const std::string& _dir, const std::string& _filePassport, const std::string& _defaults, const std::string& _switchers, const std::vector<std::string>& vars)
: dir(_dir)
{	
	std::string fileFullName = dir + _filePassport;
	std::string defaultsFileFullName = _defaults;
	std::string switchersFileFullName = _switchers;

	if (fileExistTest(fileFullName, *(defaults::defaultPinfo), *(defaults::defaultPerr), "passport") 
		&& fileExistTest(defaultsFileFullName, *(defaults::defaultPinfo), *(defaults::defaultPerr), "defaults")
		&& fileExistTest(switchersFileFullName, *(defaults::defaultPinfo), *(defaults::defaultPerr), "switchers"))
	{
		std::stringstream varsStream(StreamParser::VectorStringToString(vars));

		GetAllParamsFromParser(
			Preprocessor(fileFullName).resultStream, 
			Preprocessor(defaultsFileFullName).resultStream,
			Preprocessor(switchersFileFullName).resultStream,
			varsStream
			);

		PrintAllParams();
	}	
}//Passport(...)


//Считывание всех параметров расчета из соответствующих потоков
void Passport::GetAllParamsFromParser
(
	std::istream& mainStream, 
	std::istream& defaultStream, 
	std::istream& switcherStream, 
	std::istream& varsStream
)
{
	// 1. Разбор паспорта в целом
	
	//создаем парсер и связываем его с нужным потоком
	std::unique_ptr<StreamParser> parser;
	parser.reset(new StreamParser(mainStream, defaultStream, switcherStream, varsStream));
		
	//считываем общие парамеры
	parser->get("rho", physicalProperties.rho);
	parser->get("vInf", physicalProperties.vInf);
	parser->get("timeAccel", physicalProperties.timeAccel, &defaults::defaultTimeAccel);
	parser->get("nu", physicalProperties.nu);

	parser->get("timeStart", timeDiscretizationProperties.timeStart, &defaults::defaultTimeStart);
	parser->get("timeStop", timeDiscretizationProperties.timeStop);
	parser->get("dt", timeDiscretizationProperties.dt);
	parser->get("deltacntText", timeDiscretizationProperties.deltacntText, &defaults::defaultDeltacntText);
	parser->get("deltacntBinary", timeDiscretizationProperties.deltacntBinary, &defaults::defaultDeltacntBinary);


	parser->get("eps", wakeDiscretizationProperties.eps);
	wakeDiscretizationProperties.eps2 = wakeDiscretizationProperties.eps* wakeDiscretizationProperties.eps;
	parser->get("epscol", wakeDiscretizationProperties.epscol);
	parser->get("distKill", wakeDiscretizationProperties.distKill, &defaults::defaultDistKill);
	parser->get("delta", wakeDiscretizationProperties.delta, &defaults::defaultDelta);
	
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
		std::pair<std::string, std::string> airfoilLine = StreamParser::SplitString(airfoil[i]);

		AirfoilParams prm;
		//первая подстрока - имя файла
		prm.fileAirfoil = airfoilLine.first;

		//вторую подстроку разделяем на вектор из строк по запятым, стоящим вне фигурных скобок		
		std::vector<std::string> vecAirfoilLineSecond = StreamParser::StringToVector(airfoilLine.second, '{', '}');
				
		//создаем парсер и связываем его с параметрами профиля
		std::stringstream aflStream(StreamParser::VectorStringToString(vecAirfoilLineSecond));
		std::unique_ptr<StreamParser> parserAirfoil;

		parserAirfoil.reset(new StreamParser(aflStream, defaultStream, switcherStream, varsStream));

		//считываем нужные параметры с учетом default-значений
		parserAirfoil->get("basePoint", prm.basePoint, &defaults::defaultBasePoint);
		parserAirfoil->get("scale", prm.scale, &defaults::defaultScale);
		parserAirfoil->get("angle", prm.angle, &defaults::defaultAngle);
		parserAirfoil->get("panelsType", prm.panelsType, &defaults::defaultPanelsType);
		parserAirfoil->get("boundaryConditionSatisfaction", prm.boundaryCondition, &defaults::defaultBoundaryCondition);
		parserAirfoil->get("mechnicalSystem", prm.mechanicalSystem, &defaults::defaultMechanicalSystem);

		//отправляем считанные параметры профиля в структуру данных паспорта 
		airfoilParams.push_back(prm);

	} //for i
}//GetAllParamsFromParser(...)


//Печать всех параметров расчета в поток логов
void Passport::PrintAllParams()
{
	const std::string str = "passport info: ";
	std::ostream& out = *(defaults::defaultPinfo);

	out << str << "rho = " << physicalProperties.rho << std::endl;
	out << str << "vInf = " << physicalProperties.vInf << std::endl;
	out << str << "nu = " << physicalProperties.nu << std::endl;
	out << str << "timeStart = " << timeDiscretizationProperties.timeStart << std::endl;
	out << str << "timeStop = " << timeDiscretizationProperties.timeStop << std::endl;
	out << str << "dt = " << timeDiscretizationProperties.dt << std::endl;
	out << str << "deltacntText = " << timeDiscretizationProperties.deltacntText << std::endl;
	out << str << "deltacntBinary = " << timeDiscretizationProperties.deltacntBinary << std::endl;
	out << str << "eps = " << wakeDiscretizationProperties.eps << std::endl;
	out << str << "eps2 = " << wakeDiscretizationProperties.eps2 << std::endl;
	out << str << "epscol = " << wakeDiscretizationProperties.epscol << std::endl;
	out << str << "distKill = " << wakeDiscretizationProperties.distKill << std::endl;
	out << str << "delta = " << wakeDiscretizationProperties.delta << std::endl;
	out << str << "linearSystemSolver = " << numericalSchemes.linearSystemSolver << std::endl;
	out << str << "velocityComputation = " << numericalSchemes.velocityComputation << std::endl;
	out << str << "wakeMotionIntegrator = " << numericalSchemes.wakeMotionIntegrator << std::endl;
	
	out << str << "airfoilsDir = " << airfoilsDir << std::endl;
	out << str << "wakesDir = " << wakesDir << std::endl;

	out << str << "number of airfoils = " << airfoilParams.size() << std::endl;
	for (size_t q = 0; q < airfoilParams.size(); ++q)
	{
		out << str << "airfoil[" << q << "]_file = " << airfoilParams[q].fileAirfoil << std::endl;
		out << str << "airfoil[" << q << "]_basePoint = " << airfoilParams[q].basePoint << std::endl;
		out << str << "airfoil[" << q << "]_scale = " << airfoilParams[q].scale << std::endl;
		out << str << "airfoil[" << q << "]_angle = " << airfoilParams[q].angle << std::endl;
		out << str << "airfoil[" << q << "]_panelType = " << airfoilParams[q].panelsType << std::endl;
		out << str << "airfoil[" << q << "]_boundaryCondition = " << airfoilParams[q].boundaryCondition << std::endl;
		out << str << "airfoil[" << q << "]_mechanicalSystem = " << airfoilParams[q].mechanicalSystem << std::endl;
	}

	out << str << "fileWake = " << wakeDiscretizationProperties.fileWake << std::endl;
}//PrintAllParams()
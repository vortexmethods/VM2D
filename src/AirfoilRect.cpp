/*!
\file
\brief Файл кода с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "AirfoilRect.h"


//Считывание профиля из файла
void AirfoilRect::ReadFromFile(const std::string& dir, const AirfoilParams& param) //загрузка профиля из файла, его поворот и масштабирование
{
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;
	if (fileExistTest(filename, *(defaults::defaultPinfo), *(defaults::defaultPerr), "airfoil"))
	{
		airfoilFile.open(filename);
		StreamParser airfoilParser(airfoilFile);
		airfoilFile.close();

		rcm = param.basePoint;
		m = 0.0; //TODO
		J = 0.0; //TODO
		
		airfoilParser.get("np", np);

		airfoilParser.get("r", r);
		r.push_back(r[0]);

		for (size_t q = 0; q < r.size(); ++q)
			v.push_back({ 0.0, 0.0 });

		Move(rcm);
		Scale(param.scale);
		Rotate(param.angle);
		//в конце Rotate нормали, касательные и длины вычисляются сами
	}
}//ReadFromFile(...)

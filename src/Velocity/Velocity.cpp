/*!
\file
\brief Файл кода с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/ 


#include "Velocity.h"

//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcConvVelo()
{
	CalcConvVeloToSetOfPoints(wake.vtx, wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake);

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < epsastWake.size(); ++i)
	//	veloFile << epsastWake[i] << std::endl;
	//veloFile.close();

	//TODO пока только один профиль	
	if (boundary.size() > 0)
	{
		CalcConvVeloToSetOfPoints(boundary[0]->virtualWake, virtualVortexesParams[0].convVelo, virtualVortexesParams[0].epsastWake);
		/*std::ostringstream sss;
		sss << "virtVelo_";
		std::ofstream virtVeloFile(sss.str());
		for (int i = 0; i < convVirtualVelo[0].size(); ++i)
			virtVeloFile << convVirtualVelo[0][i] << std::endl;
		virtVeloFile.close();*/
	}
}//CalcConvVelo()


//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVelo()
{
	CalcDiffVeloToSetOfPoints(wake.vtx, wakeVortexesParams.epsastWake, wake.vtx, wakeVortexesParams.diffVelo);

	//TODO пока только один профиль	
	if (boundary.size() > 0)
	{
		CalcDiffVeloToSetOfPoints(wake.vtx, wakeVortexesParams.epsastWake, boundary[0]->virtualWake, wakeVortexesParams.diffVelo);

		CalcDiffVeloToSetOfPoints(boundary[0]->virtualWake, virtualVortexesParams[0].epsastWake, boundary[0]->virtualWake, virtualVortexesParams[0].diffVelo);

		CalcDiffVeloToSetOfPoints(boundary[0]->virtualWake, virtualVortexesParams[0].epsastWake, wake.vtx, virtualVortexesParams[0].diffVelo);
	}
}//CalcDiffVelo()

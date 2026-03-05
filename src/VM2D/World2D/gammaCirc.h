// Вычисление количества завихренности, находящейся в круге заданного (довольно большого радиуса)
//
// Следует сделать #include "gammaCirc.h" в src/VM2D/World2D/World2D.cpp 
// в функции Step() после строки, в которой вызывается CalcVortexVelo()
//
// Ниже числа в строке #define radiuses будут определять, по окружностям
// каких радиусов измеряется циркуляция поля скоростей
//
// Результат сохраняется в файл circ в рабочий каталог задачи

//#define radiuses { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
//				        2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 
//                        3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 
//                        4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0 }						
#define radiuses { 0.5, 1.0, 1.5, 2.0, 2.5 }
						

//#define scales {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5}
#define scales { 1.0 }


if (currentStep % 1 == 0)
{

	std::vector<double> rad = radiuses;
	std::vector<double> scale = scales;


	std::string circFileName = getPassport().dir + "circ";


	if (currentStep == 0)
	{
		std::ofstream circFile(circFileName.c_str());

		VMlib::PrintLogoToTextFile(circFile, circFileName, "Circulations along the circles of fixed radius");

		std::stringstream ssCirc;
		ssCirc << "step  time  scale";
		for (auto r : rad)
			ssCirc << "  R=" << r;

		VMlib::PrintHeaderToTextFile(circFile, ssCirc.str());

		circFile << std::endl;
		circFile.close();
		circFile.clear();
	}

	for (int c = 0; c < scale.size(); ++c)
	{
		std::vector<double> circRad(rad.size(), 0.0);

		int nPts = 50; //50

		double delta = 6.0; //6.0

		double R = 0.0;
		double L = 0.0;
		double epsAst = 0.0;
		double gam = 0.0;
		double curCirc = 0.0;

		for (size_t i = 0; i < rad.size(); ++i)
		{
			R = rad[i];
			curCirc = 0.0;

#pragma omp parallel for default(none) private(L,gam,epsAst) shared(R,nPts,delta,c,scale) reduction(+:curCirc) schedule(dynamic,100)
			for (int j = 0; j < (int)wake->vtx.size(); ++j)
			{
				L = wake->vtx[j].r().length();
				gam = wake->vtx[j].g();
				epsAst = scale[c] * velocity->wakeVortexesParams.epsastWake[j];

				if (R - L > delta * epsAst)
					curCirc += gam;
				else if (R - L > -delta * epsAst)
				{					
					double h = fabs(R - L);
					if (h / R > 1e-3)
					{
						double phi = atan(h / sqrt(delta * delta * epsAst * epsAst - h * h));
						double deltaPhi = (PI - 2.0 * phi) / nPts;

						double addCirc = 0.0;

						auto func = [h, epsAst, delta](double phi) -> double
							{
								return (IDPI / epsAst) * (exp(h / (epsAst * sin(phi))) * (epsAst - h / sin(phi)) - exp(-delta) * (1 + delta));
							};

						for (int s = 0; s < nPts; ++s)
						{
							double valAngle = -phi - deltaPhi * (s + 0.5);

							addCirc += func(valAngle) * deltaPhi;
						}
						addCirc *= gam;

						if (L > R)
							curCirc += addCirc;
						else
							curCirc += (gam - addCirc);
					}
					else
						curCirc += 0.5 * gam;
				}//else
			}//for j

			circRad[i] = curCirc;
		}//for i


		std::ofstream circFile(circFileName.c_str(), std::ios::app);

		circFile << currentStep << " " << getCurrentTime() << " " << scale[c];

		for (auto gamma : circRad)
			circFile << " " << gamma;
		circFile << std::endl;


		circFile.close();
		circFile.clear();

	}// for c

}
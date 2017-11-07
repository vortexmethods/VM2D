/*!
\file
\brief Файл кода с описанием класса Boundary
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Boundary.h"



//Конструктор
Boundary::Boundary(const Passport& passport_, const Airfoil& afl_, const std::vector<std::unique_ptr<Boundary>>& allBoundary_, int sheetDim_, const Wake& wake_, const Parallel& parallel_)
	: passport(passport_ ), afl(afl_), allBoundary(allBoundary_), sheetDim(sheetDim_), wake(wake_), parallel(parallel_), CC(afl.r)
{
	sheets.SetLayersDim(afl.np, sheetDim);
	virtualWake.resize(0);
}//Boundary(...)



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
Boundary::Boundary(const std::unique_ptr<Airfoil>& afl_, int sheetDim_, const Wake& wake_, const Parallel& parallel_)
	: afl(afl_), sheetDim(sheetDim_), wake(wake_), parallel(parallel_), CC(afl->r)
{
	sheets.SetLayersDim(afl->np, sheetDim);
}//Boundary(...)
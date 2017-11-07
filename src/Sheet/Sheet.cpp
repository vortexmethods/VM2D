/*!
\file
\brief Файл кода с описанием класса Sheet
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Sheet.h"


//Установка pазмерностей всех векторов и их обнуление
void Sheet::SetLayersDim(size_t np, size_t layerDim)
{
	//+1 для удобства, последний равен нулевому
	freeVortexSheet.resize(np+1);
	attachedVortexSheet.resize(np+1);
	attachedSourceSheet.resize(np+1);

	for (size_t j = 0; j < np+1; ++j)
	{
		freeVortexSheet[j].resize(layerDim, 0.0);
		attachedVortexSheet[j].resize(layerDim, 0.0);
		attachedSourceSheet[j].resize(layerDim, 0.0);
	}
};//SetLayersDim(...)
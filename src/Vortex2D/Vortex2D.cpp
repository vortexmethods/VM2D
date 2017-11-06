/*!
\file
\brief Файл кода с описанием класса Vortex2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "Vortex2D.h"


MPI_Datatype Vortex2D::mpiVortex2D;


void Vortex2D::CreateMpiType()
{
	int          len[4] = { 1, 2, 1, 1 };
	MPI_Aint     pos[4] = { 0, offsetof(Vortex2D, pos), offsetof(Vortex2D, gam), sizeof(Vortex2D) };
	MPI_Datatype typ[4] = { MPI_LB, MPI_DOUBLE, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(4, len, pos, typ, &mpiVortex2D);
	MPI_Type_commit(&mpiVortex2D);
}//CreateMpiType()
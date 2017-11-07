#include "MechanicsRigidImmovable.h"

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidImmovable::GetHydroDynamForce()
{
	hydroDynamForce = { 0.0, 0.0 };

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k

	for (size_t i = 0; i < afl.np; ++i)
	{
		double GamK = boundary.virtualWake[i].g();
		double deltaK =  boundary.sheets.freeVortexSheet[i][0] * afl.len[i] - afl.gammaThrough[i];  		 //afl.gammaThrough[i];
		Point2D VelK = virtVortParams.convVelo[i] /* + virtVortParams.diffVelo[i]*/ + passport.physicalProperties.V0();
		Point2D rK =  0.5 * (afl.r[i + 1] + afl.r[i]);	//boundary.virtualWake[i].r();		

		hDFGam += GamK * Point2D({ VelK[1], -VelK[0] });
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
	}

	hydroDynamForce = hDFGam + (1.0 / passport.timeDiscretizationProperties.dt) * hDFdelta;
}
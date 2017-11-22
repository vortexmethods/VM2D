/*!
\file
\brief В РАБОТЕ Файл кода с описанием класса VelocityFourier
\warning Пока в состоянии нерабочего прототипа
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "VelocityFourier.h"

 // TODO проверить resize

//пересчет глобального номера ячейки по "локальным координатам" ячейки
inline int VelocityFourier::cellIj(int i, int j) const
{ 
	return j*nCell[0] + i; 
}

//"размазанная" завихренность
//inline double VelocityFourier::W(const Point2D& pt) const
//{
//	return M4(pt[0] / h[0]) * M4(pt[1] / h[1]);
//}

//заполнение списков в массиве vic
void VelocityFourier::FillVic()
{
	vic.resize(nCell[0] * nCell[1]);
	for (size_t q = 0; q < wake.vtx.size(); ++q)
	{
		const Point2D& pos = wake.vtx[q].r();
		int i = static_cast<int>((pos[0] - corner[0]) / h[0]);
		int j = static_cast<int>((pos[1] - corner[1]) / h[1]);

		vic[cellIj(i, j)].push_back(q);
	}
}

//завихренность, приведенная к узлу (i,j)
double VelocityFourier::w(int i, int j) const
{
	int leftCell = std::max(0, i - 2);
	int rightCell = std::min(nCell[0], i + 2);

	int downCell = std::max(0, j - 2);
	int upCell = std::min(nCell[1], j + 2);

	double res = 0.0;
	Point2D posNode = corner + Point2D({ i*h[0], j*h[1] });

	for (int i = leftCell; i < rightCell; ++i)
	for (int j = downCell; j < upCell; ++j)
	{

		int cell = cellIj(i, j);
		for (size_t q = 0; q < vic[cell].size(); ++q)
		{
			const int& globQ = vic[cell][q];
			const Vortex2D& v = wake.vtx[globQ];
			
			res += Kernel((v.r()[0] - posNode[0]) / h[0]) * Kernel((v.r()[1] - posNode[1]) / h[1]);
		}
	}
	return res / (h[0] * h[1]);
}


/*
void VelocityFourier::CalcConvVelo(double dt)
{
	/// \todo Доделать OMP, MPI
	const int& np = parallel.nProcWork;
	const int& id = parallel.myidWork;

	int nvPerP = wake.nv / np;

	std::vector<int> len;
	std::vector<int> disp;

	FillVic();
	Eigen::MatrixXd omega, bkx, bky;
	omega.setZero(2 * nNode[0], 2 * nNode[1]);
	bkx.setZero(2 * nNode[0], 2 * nNode[1]);
	bky.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int i = 0; i < nNode[0]; ++i)
		for (int j = 0; j < nNode[1]; ++j)
		{
			omega(i, j) = w(i, j);
			Point2D skos = Skos(corner + Point2D({ i*h[0], j*h[1] }), { 0.0, 0.0 });
			bkx(i, j) = skos[0];
			bky(i, j) = skos[1];
		}

	for (int i = 0; i < nNode[0] - 1; i++)
		for (int j = 0; j < nNode[1]; j++)
		{
			bkx(2 * nNode[0] - 1 - i, j) = bkx(i + 1, j);
			bky(2 * nNode[0] - 1 - i, j) = -bky(i + 1, j);
		}

	for (int i = 0; i < nNode[0]; i++)
		for (int j = 0; j < nNode[1] - 1; j++)
		{
			bkx(i, 2 * nNode[1] - 1 - j) = -bkx(i, j + 1);
			bky(i, 2 * nNode[1] - 1 - j) = bky(i, j + 1);
		}

	for (int i = 0; i < nNode[0] - 1; i++)
		for (int j = 0; j < nNode[1] - 1; j++)
		{
			bkx(2 * nNode[0] - 1 - i, 2 * nNode[1] - 1 - j) = -bkx(i + 1, j + 1);
			bky(2 * nNode[0] - 1 - i, 2 * nNode[1] - 1 - j) = -bky(i + 1, j + 1);
		}

	if (id == 0)
		convVelo.resize(wake.nv);

	std::ofstream strOmega("matrOmega");
	SaveToStream(omega, strOmega);
	strOmega.close();

	std::ofstream strB("matrB");
	SaveToStream(bkx, strB);
	strB.close();


	//FFT
	Eigen::FFT<double> fft;

	Eigen::MatrixXcd fOmega;
	fOmega.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int k = 0; k < omega.rows(); k++) {
		Eigen::VectorXcd tmpOmega(2 * nNode[0]);
		fft.fwd(tmpOmega, omega.row(k));
		fOmega.row(k) = tmpOmega;
	}

	for (int k = 0; k < omega.cols(); k++) {
		Eigen::VectorXcd tmpOmega(2 * nNode[1]);
		fft.fwd(tmpOmega, fOmega.col(k));
		fOmega.col(k) = tmpOmega;
	}

	Eigen::MatrixXcd fBkx;
	fBkx.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int k = 0; k < bkx.rows(); k++) {
		Eigen::VectorXcd tmpBkx(2 * nNode[0]);
		fft.fwd(tmpBkx, bkx.row(k));
		fBkx.row(k) = tmpBkx;
	}

	for (int k = 0; k < bkx.cols(); k++) {
		Eigen::VectorXcd tmpBkx(2 * nNode[1]);
		fft.fwd(tmpBkx, fBkx.col(k));
		fBkx.col(k) = tmpBkx;
	}

	Eigen::MatrixXcd fBky;
	fBky.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int k = 0; k < bky.rows(); k++) {
		Eigen::VectorXcd tmpBky(2 * nNode[0]);
		fft.fwd(tmpBky, bky.row(k));
		fBky.row(k) = tmpBky;
	}

	for (int k = 0; k < bky.cols(); k++) {
		Eigen::VectorXcd tmpBky(2 * nNode[1]);
		fft.fwd(tmpBky, fBky.col(k));
		fBky.col(k) = tmpBky;
	}

	std::ofstream strFOmega("matrFOmega");
	SaveToStream(fOmega, strFOmega);
	strFOmega.close();

	std::ofstream strFB("matrFB");
	SaveToStream(fBky, strFB);
	strFB.close();

	Eigen::MatrixXcd fUx, fUy;
	fUx.setZero(2 * nNode[0], 2 * nNode[1]);
	fUy.setZero(2 * nNode[0], 2 * nNode[1]);
	for (int i = 0; i < 2 * nNode[0]; ++i)
		for (int j = 0; j < 2 * nNode[1]; ++j)
		{
			fUx(i, j) = fBkx(i, j)*fOmega(i, j);
			fUy(i, j) = fBky(i, j)*fOmega(i, j);
		}

	Eigen::MatrixXcd ux;
	ux.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int k = 0; k < fUx.rows(); k++) {
		Eigen::VectorXcd tmpUx(2 * nNode[0]);
		fft.inv(tmpUx, fUx.row(k));
		ux.row(k) = tmpUx;
	}

	for (int k = 0; k < fUx.cols(); k++) {
		Eigen::VectorXcd tmpUx(2 * nNode[1]);
		fft.inv(tmpUx, ux.col(k));
		ux.col(k) = tmpUx;
	}

	Eigen::MatrixXcd uy;
	uy.setZero(2 * nNode[0], 2 * nNode[1]);

	for (int k = 0; k < fUy.rows(); k++) {
		Eigen::VectorXcd tmpUy(2 * nNode[0]);
		fft.inv(tmpUy, fUy.row(k));
		uy.row(k) = tmpUy;
	}

	for (int k = 0; k < fUy.cols(); k++) {
		Eigen::VectorXcd tmpUy(2 * nNode[1]);
		fft.inv(tmpUy, uy.col(k));
		uy.col(k) = tmpUy;
	}

	std::vector<Point2D> uxy;
	for (int i = 0; i < nNode[0]; ++i)
		for (int j = 0; j < nNode[1]; ++j)
			uxy.push_back(1.0 / (sqrt(nNode[0] * nNode[1] * (nNode[0] - 2) * (nNode[1] - 2))) * Point2D({ ux(i, j).real(), uy(i, j).real() }));

	std::ofstream strU("matrU");
	SaveToStream(uxy, strU);
	strU.close();

	//MPI_Gatherv(locConvVelo.data(), len[id], Point2D::mpiPoint2D, convVelo.data(), len.data(), disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
}


*/
#pragma once

#include <vector>

class TriDiagMatrix
{
public:
	std::vector<double> a, b, c;
	bool cyclic;

	TriDiagMatrix(size_t dim)
	{
		a.resize(dim);
		b.resize(dim);
		c.resize(dim);
		cyclic = false;
	};

	~TriDiagMatrix() {};

	void setABC(const std::vector<double>& a_, const std::vector<double>& b_, const std::vector<double>& c_)
	{
		a = a_;
		b = b_;
		c = c_;
	};

	double& operator()(size_t i, size_t j)
	{
#pragma warning (push)
#pragma warning (disable: 4715)
		if (i == j) return b[i];
		else if (i == ((j + 1) % a.size())) return a[i];
		else if (((i + 1) % a.size()) == j) return c[i];
		else return b[i];
#pragma warning (pop)
	}

	std::vector<double> solve(const std::vector<double>& d) const
	{
		if (cyclic)
			return cyclicTriDiagRun(d);
		else
			return TriDiagRun(d);
	}

	std::vector<double> TriDiagRun(const std::vector<double>& d) const
	{
		const size_t dim = a.size();

		std::vector<double> x(dim, 0.0);

		std::vector<double> alpha(dim - 1, 0.0), beta(dim - 1, 0.0);

		alpha[0] = -c[0] / b[0];
		beta[0] = d[0] / b[0];

		double den;

		for (size_t i = 1; i + 1 < dim; ++i)
		{
			den = b[i] + a[i] * alpha[i - 1];

			alpha[i] = -c[i] / den;
			beta[i] = (d[i] - a[i] * beta[i - 1]) / den;
		}

		x[dim - 1] = (d[dim - 1] - a[dim - 1] * beta[dim - 2]) / (b[dim - 1] + a[dim - 1] * alpha[dim - 2]);

		for (size_t i = dim - 2; i != -1; --i)
			x[i] = alpha[i] * x[i + 1] + beta[i];

		return x;
	}

	std::vector<double> cyclicTriDiagRun(const std::vector<double>& d) const
	{
		const size_t dim = a.size();

		std::vector<double> x(dim, 0.0);

		std::vector<double> alpha(dim - 1, 0.0), beta(dim - 1, 0.0), gamma(dim - 1, 0.0);

		alpha[0] = -c[0] / b[0];
		beta[0] = -a[0] / b[0];
		gamma[0] = d[0] / b[0];

		double den;

		for (size_t i = 1; i + 1 < dim; ++i)
		{
			den = b[i] + a[i] * alpha[i - 1];

			alpha[i] = -c[i] / den;
			beta[i] = -a[i] * beta[i - 1] / den;
			gamma[i] = (d[i] - a[i] * gamma[i - 1]) / den;
		}

		std::vector<double> lambda(dim - 1, 0.0), mu(dim - 1, 0.0);


		lambda[dim - 2] = alpha[dim - 2] + beta[dim - 2];
		mu[dim - 2] = gamma[dim - 2];

		for (size_t i = dim - 3; i != -1; --i)
		{
			lambda[i] = alpha[i] * lambda[i + 1] + beta[i];
			mu[i] = alpha[i] * mu[i + 1] + gamma[i];
		}

		x[dim - 1] = (d[dim - 1] - c[dim - 1] * mu[0] - a[dim - 1] * gamma[dim - 2]) / (b[dim - 1] + c[dim - 1] * lambda[0] + a[dim - 1] * (alpha[dim - 2] + beta[dim - 2]));
		for (size_t i = 0; i + 1 < dim; ++i)
		{
			x[i] = lambda[i] * x[dim - 1] + mu[i];
		}

		return x;
	}
};
//#pragma once

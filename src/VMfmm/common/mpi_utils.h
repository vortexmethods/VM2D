#pragma once
#include "../common/utils.h"
#include <utility>
#include <iostream>

#ifndef __NVCC__
#include <oneapi/tbb/parallel_for.h>
#endif

#ifdef FMM_MPI
#include "mpi.h"

namespace fmm {

inline int RootID()
{
	return 0;
}

inline int MyID()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank;
}

inline bool IAmRoot()
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank == RootID();
}

inline int NProc()
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	return size;
}

inline std::pair<size_t, size_t> LocalPart(size_t begin, size_t end)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int N = NProc();
	size_t size = end - begin;
	return { rank * size / N, (rank + 1) * size / N };
}

inline int LocalPart(size_t begin, size_t end, int rank)
{
	int N = NProc();
	size_t size = end - begin;
	return (rank + 1) * size / N - rank * size / N;
}

#ifndef __NVCC__
template <typename T>
inline void AllReduce(T* data, int size)
{
	MPI_Datatype mpi_type;
	if constexpr (std::is_same_v<T, double>)
		mpi_type = MPI_DOUBLE;
	if constexpr (std::is_same_v<T, Vector3d>)
	{
		MPI_Type_contiguous(3, MPI_DOUBLE, &mpi_type);
		MPI_Type_commit(&mpi_type);
	}
	if constexpr (std::is_same_v<T, std::complex<double>>)
		mpi_type = MPI_COMPLEX16;

	const int nproc = NProc();
	std::vector<int> send_displs(nproc + 1);
	for (int i = 0; i <= nproc; ++i)
		send_displs[i] = size * i / nproc;
	std::vector<int> send_counts(nproc);
	for (int i = 0; i < nproc; ++i)
		send_counts[i] = send_displs[i + 1] - send_displs[i];

	const int local_size = send_counts[MyID()];
	std::vector<int> recv_displs(nproc, 0);
	for (int i = 1; i < nproc; ++i)
		recv_displs[i] = recv_displs[i - 1] + local_size;
	std::vector<int> recv_counts(nproc, local_size);

	std::unique_ptr<T[]> buf(new T[local_size * nproc]);
	MPI_Alltoallv(data, send_counts.data(), send_displs.data(), mpi_type, buf.get(), recv_counts.data(), recv_displs.data(), mpi_type, MPI_COMM_WORLD);
	tbb::parallel_for(0, local_size, [&](int i) {
		for (int j = 1; j < nproc; ++j)
			buf[i] += buf[i + j * local_size];
		});
	MPI_Allgatherv(buf.get(), local_size, mpi_type, data, send_counts.data(), send_displs.data(), mpi_type, MPI_COMM_WORLD);
}
#endif

}
#else

namespace fmm {

inline bool IAmRoot()
{
	return true;
}

}

#endif

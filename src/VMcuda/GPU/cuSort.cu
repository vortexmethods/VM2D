/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.5    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/08/29     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: cuSort.cu                                                        |
| Info: Source code of BHgpu                                                  |
|                                                                             |
| This file is part of BHgpu.                                                 |
| BHcu is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BHgpu.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Реализация CUDA-ядер для сортировки на GPU
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#include "cuSort.cuh"

//#include "cub/cub.cuh"   // or equivalently 
#include <cub/device/device_radix_sort.cuh>


namespace BHcu
{

    void RadixSortFromCUB(        
        const int* MmortonCodesKeyUnsortd, int* MmortonCodesKeyd, \
        const int* MmortonCodesIdxUnsortd, int* MmortonCodesIdxd,
        int num_items, int begin_bit, int end_bit)
    {
        ///RadixSort

            // Declare, allocate, and initialize device-accessible pointers for sorting data
            //int  num_items = 7;      // e.g., 7
            //int  *d_keys_in;         // e.g., [8, 6, 7, 5, 3, 0, 9]
            //int  *d_keys_out;        // e.g., [        ...        ]
            //int  *d_values_in;       // e.g., [0, 1, 2, 3, 4, 5, 6]
            //int  *d_values_out;      // e.g., [        ...        ]

            //d_keys_in = (int*)cudaNew(7, sizeof(int));
            //d_keys_out = (int*)cudaNew(7, sizeof(int));
            //d_values_in = (int*)cudaNew(7, sizeof(int));
            //d_values_out = (int*)cudaNew(7, sizeof(int));

            //int h_keys_in[] = { 8, 6, 7, 5, 3, 0, 9 };
            //int h_values_in[] = { 0, 1, 2, 3, 4, 5, 6 };
            //int h_keys_out[7] ;
            //int h_values_out[7];

            //cudaCopyVecToDevice(h_keys_in, d_keys_in, 7, sizeof(int));
            //cudaCopyVecToDevice(h_values_in, d_values_in, 7, sizeof(int));


            // Determine temporary device storage requirements
        void* d_temp_storage = NULL;
        size_t   temp_storage_bytes = 0;
        //cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);

        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, \
            MmortonCodesKeyUnsortd, MmortonCodesKeyd, \
            MmortonCodesIdxUnsortd, MmortonCodesIdxd, \
            num_items, begin_bit, end_bit);

        // Allocate temporary storage
        cudaMalloc(&d_temp_storage, temp_storage_bytes);

        // Run sorting operation
        //cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, d_keys_in, d_keys_out, d_values_in, d_values_out, num_items);
        cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes, \
            MmortonCodesKeyUnsortd, MmortonCodesKeyd, \
            MmortonCodesIdxUnsortd, MmortonCodesIdxd, \
            num_items, begin_bit, end_bit);

        // d_keys_out            <-- [0, 3, 5, 6, 7, 8, 9]
        // d_values_out          <-- [5, 4, 3, 1, 2, 0, 6]

        //cudaCopyVecFromDevice(d_keys_out, h_keys_out, 7, sizeof(int));
        //cudaCopyVecFromDevice(d_values_out, h_values_out, 7, sizeof(int));

        //for (int i = 0; i < num_items; ++i)
        //	std::cout << h_keys_out[i] << " ";
        //std::cout << "\n";
        //for (int i = 0; i < num_items; ++i)
        //	std::cout << h_values_out[i] << " ";

        cudaFree(d_temp_storage);

    }

}
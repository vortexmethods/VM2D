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
| File name: operations.cuh                                                   |
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
\brief Вспомогательные операции
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/


#ifndef MLP_H_
#define MLP_H_

#include <cuda.h>
#include <cudnn.h>


enum class ConvAlgoSelector {
    AUTO,
    IMPLICIT_GEMM,
    IMPLICIT_PRECOMP_GEMM,
    GEMM,
    DIRECT,
    FFT,
    FFT_TILING,
    WINOGRAD,
    WINOGRAD_NONFUSED
};

struct Layer {
    int input_size;
    int output_size;

    float* d_weights = nullptr;
    float* d_biases = nullptr;
    float* d_output = nullptr;

    cudnnTensorDescriptor_t input_desc;
    cudnnTensorDescriptor_t output_desc;
    cudnnTensorDescriptor_t bias_desc;
    cudnnFilterDescriptor_t weight_desc;
    cudnnActivationDescriptor_t activation_desc;
    cudnnConvolutionDescriptor_t conv_desc;

    cudnnConvolutionFwdAlgo_t conv_algo;
    void* d_workspace = nullptr;
    size_t workspace_size = 0;

    bool has_activation = false; // 👈 добавили флаг

    void destroy() {
        cudnnDestroyTensorDescriptor(input_desc);
        cudnnDestroyTensorDescriptor(output_desc);
        cudnnDestroyTensorDescriptor(bias_desc);
        cudnnDestroyFilterDescriptor(weight_desc);
        cudnnDestroyConvolutionDescriptor(conv_desc);
        if (has_activation) {
            cudnnDestroyActivationDescriptor(activation_desc);
        }
        cudaFree(d_workspace);
        cudaFree(d_weights);
        cudaFree(d_biases);
        cudaFree(d_output);
    }

};



class CudaNN 
{
public:
	std::vector<std::vector<float>> weights;
	std::vector<std::vector<float>> biases;
	std::vector<int> input_size;
	std::vector<int> output_size;

    cudnnHandle_t cudnn;
    std::vector<Layer> network;

	CudaNN();	
	~CudaNN() = default;

    void read_params(const std::string& filename, int layer, int input_size, int output_size);

    int start(int batch_size, float* h_input, float* h_output);

	void setup_layer(Layer& layer, const std::string& layer_name,
		int input_size, int output_size, int batch_size,
		ConvAlgoSelector algo_choice = ConvAlgoSelector::AUTO,
		bool with_activation = true);


	std::string make_layer_signature(const std::string& layer_name,
		int batch_size,
		int input_size,
		int output_size,
		cudnnConvolutionFwdAlgo_t algo,
		size_t workspace_size);

    cudnnConvolutionFwdAlgo_t resolve_algo(ConvAlgoSelector sel) {
        switch (sel) {
        case ConvAlgoSelector::IMPLICIT_GEMM: return CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_GEMM;
        case ConvAlgoSelector::IMPLICIT_PRECOMP_GEMM: return CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_PRECOMP_GEMM;
        case ConvAlgoSelector::GEMM: return CUDNN_CONVOLUTION_FWD_ALGO_GEMM;
        case ConvAlgoSelector::DIRECT: return CUDNN_CONVOLUTION_FWD_ALGO_DIRECT;
        case ConvAlgoSelector::FFT: return CUDNN_CONVOLUTION_FWD_ALGO_FFT;
        case ConvAlgoSelector::FFT_TILING: return CUDNN_CONVOLUTION_FWD_ALGO_FFT_TILING;
        case ConvAlgoSelector::WINOGRAD: return CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD;
        case ConvAlgoSelector::WINOGRAD_NONFUSED: return CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD_NONFUSED;
        default: return CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_GEMM;  // fallback
        }
    };

    void setup_layers(int batch_size);
    void run_layer(Layer& layer, float* d_input, const char* layer_name);
    void destroy_layers();    
};

//#include <cuda.h>

//int cudnnStart(int batch_size, float* h_input, float* h_output);



#endif
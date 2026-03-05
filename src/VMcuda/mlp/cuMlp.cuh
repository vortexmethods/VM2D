/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cuMlp.cuh                                                        |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовки интерфейса для работы с многослойным перцептроном на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Сухова Екатерина Юрьевна
\Version 1.14
\date 6 марта 2026 г.
*/


#ifndef MLP_H
#define MLP_H

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
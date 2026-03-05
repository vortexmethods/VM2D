/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cuMlp.cu                                                         |
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
\brief Реализация функций для работы с многослойным перцептроном на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Сухова Екатерина Юрьевна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifdef USE_CUDNN

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cuda_runtime.h>
#include <cudnn.h>

#include "mlp_cudnn.cuh"

#define checkCUDNN(expr) do { \
    cudnnStatus_t status = (expr); \
    if (status != CUDNN_STATUS_SUCCESS) { \
        std::cerr << "cuDNN error: " << cudnnGetErrorString(status) << std::endl; \
        std::exit(EXIT_FAILURE); \
    } \
} while (0)

const char* cudnnGetConvolutionFwdAlgoString(cudnnConvolutionFwdAlgo_t algo) {
    switch (algo) {
    case CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_GEMM: return "IMPLICIT_GEMM";
    case CUDNN_CONVOLUTION_FWD_ALGO_IMPLICIT_PRECOMP_GEMM: return "IMPLICIT_PRECOMP_GEMM";
    case CUDNN_CONVOLUTION_FWD_ALGO_GEMM: return "GEMM";
    case CUDNN_CONVOLUTION_FWD_ALGO_DIRECT: return "DIRECT";
    case CUDNN_CONVOLUTION_FWD_ALGO_FFT: return "FFT";
    case CUDNN_CONVOLUTION_FWD_ALGO_FFT_TILING: return "FFT_TILING";
    case CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD: return "WINOGRAD";
    case CUDNN_CONVOLUTION_FWD_ALGO_WINOGRAD_NONFUSED: return "WINOGRAD_NONFUSED";
    case CUDNN_CONVOLUTION_FWD_ALGO_COUNT: return "COUNT";
    default: return "UNKNOWN";
    }
}



void CudaNN::read_params(const std::string & filename, int layer, int input_size, int output_size) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Cannot open weights file: " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    weights[layer].resize(input_size * output_size);
    biases[layer].resize(output_size);
    file.read(reinterpret_cast<char*>(weights[layer].data()), weights[layer].size() * sizeof(float));
    file.read(reinterpret_cast<char*>(biases[layer].data()), biases[layer].size() * sizeof(float));
    file.close();
}

CudaNN::CudaNN()
{
    weights.resize(3);
    biases.resize(3);
    input_size.resize(3);
    output_size.resize(3);

    read_params("../src/VMcuda/mlp/w1.bin", 0, 3, 32);
    read_params("../src/VMcuda/mlp/w2.bin", 1, 32, 32);
    read_params("../src/VMcuda/mlp/w3.bin", 2, 32, 2);

    network.resize(3);
        
    checkCUDNN(cudnnCreate(&cudnn));
}

void CudaNN::setup_layers(int batch_size)
{
    // Слой 0: вход → скрытый
    setup_layer(network[0], "input_hidden", 3, 32, batch_size, ConvAlgoSelector::IMPLICIT_PRECOMP_GEMM, true);

    // Слой 1: скрытый → скрытый
    setup_layer(network[1], "hidden_hidden", 32, 32, batch_size, ConvAlgoSelector::IMPLICIT_PRECOMP_GEMM, true);

    // Слой 2: скрытый → выход, без активации
    setup_layer(network[2], "hidden_output", 32, 2, batch_size, ConvAlgoSelector::IMPLICIT_PRECOMP_GEMM, false);      // жёсткий выбор

    cudaMemcpy(network[0].d_weights, weights[0].data(), weights[0].size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(network[0].d_biases, biases[0].data(), biases[0].size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(network[1].d_weights, weights[1].data(), weights[1].size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(network[1].d_biases, biases[1].data(), biases[1].size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(network[2].d_weights, weights[2].data(), weights[2].size() * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(network[2].d_biases, biases[2].data(), biases[2].size() * sizeof(float), cudaMemcpyHostToDevice);
}

void CudaNN::setup_layer(Layer& layer, const std::string& layer_name,
    int input_size, int output_size, int batch_size,
    ConvAlgoSelector algo_choice,
    bool with_activation)  // параметр для активации) 
{
    layer.input_size = input_size;
    layer.output_size = output_size;
    layer.has_activation = with_activation;

    // Создание дескрипторов
    checkCUDNN(cudnnCreateTensorDescriptor(&layer.input_desc));
    checkCUDNN(cudnnCreateTensorDescriptor(&layer.output_desc));
    checkCUDNN(cudnnCreateTensorDescriptor(&layer.bias_desc));
    checkCUDNN(cudnnCreateFilterDescriptor(&layer.weight_desc));
    checkCUDNN(cudnnCreateConvolutionDescriptor(&layer.conv_desc));

    // Настройка input/output
    checkCUDNN(cudnnSetTensor4dDescriptor(layer.input_desc, CUDNN_TENSOR_NCHW, CUDNN_DATA_FLOAT, batch_size, input_size, 1, 1));
    checkCUDNN(cudnnSetTensor4dDescriptor(layer.output_desc, CUDNN_TENSOR_NCHW, CUDNN_DATA_FLOAT, batch_size, output_size, 1, 1));
    checkCUDNN(cudnnSetTensor4dDescriptor(layer.bias_desc, CUDNN_TENSOR_NCHW, CUDNN_DATA_FLOAT, 1, output_size, 1, 1));
    checkCUDNN(cudnnSetFilter4dDescriptor(layer.weight_desc, CUDNN_DATA_FLOAT, CUDNN_TENSOR_NCHW, output_size, input_size, 1, 1));

    if (with_activation)
    {
        checkCUDNN(cudnnCreateActivationDescriptor(&layer.activation_desc));
        checkCUDNN(cudnnSetActivationDescriptor(layer.activation_desc, CUDNN_ACTIVATION_TANH, CUDNN_PROPAGATE_NAN, 0.0));
    }
    else
        layer.activation_desc = nullptr;

    checkCUDNN(cudnnSetConvolution2dDescriptor(layer.conv_desc, 0, 0, 1, 1, 1, 1, CUDNN_CROSS_CORRELATION, CUDNN_DATA_FLOAT));

    if (algo_choice == ConvAlgoSelector::AUTO) {
        cudnnConvolutionFwdAlgoPerf_t perf_results[10];
        int algo_count = 0;
        checkCUDNN(cudnnGetConvolutionForwardAlgorithm_v7(
            cudnn,
            layer.input_desc,
            layer.weight_desc,
            layer.conv_desc,
            layer.output_desc,
            10,
            &algo_count,
            perf_results));

        //std::cout << "[Layer setup] AUTO mode: found " << algo_count << " algorithms\n";
        //for (int i = 0; i < algo_count; ++i) {
        //    std::cout << "  Algo " << i << ": " << cudnnGetConvolutionFwdAlgoString(perf_results[i].algo)
        //        << " | Time = " << perf_results[i].time << " ms | WS = "
        //        << perf_results[i].memory / 1024.0 << " KB\n";
        //}

        layer.conv_algo = perf_results[0].algo;
    }
    else {
        layer.conv_algo = resolve_algo(algo_choice);
        //std::cout << "[Layer setup] Manual mode: using algorithm "
        //    << cudnnGetConvolutionFwdAlgoString(layer.conv_algo) << "\n";
    }

    // Workspace
    checkCUDNN(cudnnGetConvolutionForwardWorkspaceSize(
        cudnn,
        layer.input_desc,
        layer.weight_desc,
        layer.conv_desc,
        layer.output_desc,
        layer.conv_algo,
        &layer.workspace_size));

    //std::cout << "[Layer setup] Workspace size = " << layer.workspace_size / 1024.0 << " KB\n";

    cudaMalloc(&layer.d_workspace, layer.workspace_size);
    cudaMalloc(&layer.d_weights, input_size * output_size * sizeof(float));
    cudaMalloc(&layer.d_biases, output_size * sizeof(float));
    cudaMalloc(&layer.d_output, batch_size * output_size * sizeof(float));

    //std::string signature = make_layer_signature(layer_name, batch_size, input_size, output_size, layer.conv_algo, layer.workspace_size);
    //std::cout << signature << std::endl;
}

std::string CudaNN:: make_layer_signature(const std::string& layer_name,
    int batch_size,
    int input_size,
    int output_size,
    cudnnConvolutionFwdAlgo_t algo,
    size_t workspace_size) {
    std::ostringstream oss;
    oss << "[Signature] Layer: " << layer_name
        << " | Algo: " << cudnnGetConvolutionFwdAlgoString(algo)
        << " | B=" << batch_size
        << ", In=" << input_size
        << ", Out=" << output_size
        << ", WS=" << workspace_size / 1024.0 << " KB";
    return oss.str();
}



void CudaNN::run_layer(Layer& layer, float* d_input, const char* layer_name) {
    float alpha = 1.0f, beta = 0.0f;

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    checkCUDNN(cudnnConvolutionForward(
        cudnn,
        &alpha,
        layer.input_desc, d_input,
        layer.weight_desc, layer.d_weights,
        layer.conv_desc,
        layer.conv_algo,
        layer.d_workspace, layer.workspace_size,
        &beta,
        layer.output_desc, layer.d_output));

    checkCUDNN(cudnnAddTensor(
        cudnn,
        &alpha,
        layer.bias_desc, layer.d_biases,
        &alpha,
        layer.output_desc, layer.d_output));

    if (layer.has_activation) {
        checkCUDNN(cudnnActivationForward(
            cudnn,
            layer.activation_desc,
            &alpha,
            layer.output_desc, layer.d_output,
            &beta,
            layer.output_desc, layer.d_output));
    }

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float ms = 0;
    cudaEventElapsedTime(&ms, start, stop);
    //std::cout << "Layer " << layer_name << " took " << ms << " ms" << std::endl;

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

int CudaNN::start(int batch_size, float* h_input, float* h_output)
{    
    // Входной буфер (batch_size x 3)
    //float* h_input = new float[batch_size * 3];
    //for (int i = 0; i < batch_size * 3; ++i)
    //    h_input[i] = float(i % 10) / 10;

    //h_input[0] = 0.2f;
    //h_input[1] = 0.3f;
    //h_input[2] = 0.5f;

    std::vector<bool> flagY(batch_size, false);
    

    for (int i = 0; i < batch_size; ++i)
    {
        if (h_input[3 * i + 0] < 0)
            h_input[3 * i + 0] = -h_input[3 * i + 0];

        if (h_input[3 * i + 0] > 3.0f)
            h_input[3 * i + 0] = 3.0f;


        if (h_input[3 * i + 1] < 0)
        {
            h_input[3 * i + 1] = -h_input[3 * i + 1];
            flagY[i] = true;
        }
    

        if (h_input[3 * i + 1] > 3.0f)
            h_input[3 * i + 1] = 3.0f;


        if (h_input[3 * i + 2] > 1.0f)
            h_input[3 * i + 2] = 1.0f;
    }

    float* d_input;
    cudaMalloc(&d_input, batch_size * 3 * sizeof(float));
    cudaMemcpy(d_input, h_input, batch_size * 3 * sizeof(float), cudaMemcpyHostToDevice);

    run_layer(network[0], d_input, "1");
    run_layer(network[1], network[0].d_output, "2");
    run_layer(network[2], network[1].d_output, "3");

    cudaMemcpy(h_output, network[2].d_output, batch_size * 2 * sizeof(float), cudaMemcpyDeviceToHost);

    for (int i = 0; i < batch_size; ++i)
        if (flagY[i] < 0)
            h_output[2 * i + 1] = -h_output[2 * i + 1];    

    //std::cout << "Sample outputs:\n";
    //for (int i = 0; i < std::min(5, batch_size); ++i) {
    //    std::cout << "Sample[" << i << "]: "
    //              << h_output[i * 2] << ", "
    //              << h_output[i * 2 + 1] << std::endl;
    //}

    //std::cout << "Sample[" << 0 << "]: "
    //              << h_output[0 * 2] << ", "
    //              << h_output[0 * 2 + 1] << std::endl;

    cudaFree(d_input);

    return 0;
}

void CudaNN::destroy_layers() 
{
    for (auto& layer : network)
        layer.destroy();
}

#endif
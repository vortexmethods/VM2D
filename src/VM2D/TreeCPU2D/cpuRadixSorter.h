/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cpuRadixSorter.h                                                 |
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
\brief Заголовок класса для реализации поразрядной сортировки (RadixSort) на CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef CPURADIXSORT_H
#define CPURADIXSORT_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <omp.h>
#include <string>

namespace VM2D
{

    // SFINAE для проверки тривиальной копируемости (оптимизация)
    template<typename T>
    using enable_if_trivially_copyable = typename std::enable_if<std::is_trivially_copyable<T>::value>::type;

    template<typename T>
    using enable_if_not_trivially_copyable = typename std::enable_if<!std::is_trivially_copyable<T>::value>::type;

    // Основной шаблонный класс
    template<typename ValueType>
    class OptimizedRadixSorter {
    private:
        // Переиспользуемые буферы
        std::vector<int> keyBuffer;
        std::vector<ValueType> valueBuffer;

        std::vector<int> globalHist;
        std::vector<std::vector<int>> histThreads;
        std::vector<std::vector<int>> localStarts;
        std::vector<std::vector<int>> localCounters;

        int numThreads;
        size_t currentCapacity;
        bool initialized;

        constexpr static int RADIX = 256;
        constexpr static int MASK = 0xFF;
        constexpr static int PASSES = 4; // для 32-bit int

        void initializeThreadData(int threads) {
            numThreads = threads;

            // Выделяем память для гистограмм потоков
            histThreads.resize(numThreads);
            for (int t = 0; t < numThreads; ++t) {
                histThreads[t].assign(RADIX, 0);
            }

            // Выделяем память для локальных стартовых позиций
            localStarts.resize(numThreads);
            localCounters.resize(numThreads);
            for (int t = 0; t < numThreads; ++t) {
                localStarts[t].assign(RADIX, 0);
                localCounters[t].assign(RADIX, 0);
            }

            // Глобальная гистограмма
            globalHist.assign(RADIX, 0);

            initialized = true;
        }

        void ensureCapacity(size_t n) {
            if (n > currentCapacity) {
                // Увеличиваем буферы с запасом 10% для будущих сортировок
                keyBuffer.resize(n + n / 10);
                valueBuffer.resize(n + n / 10);
                currentCapacity = keyBuffer.size();
            }
        }

        void resetHistograms() {
            // Быстрый сброс гистограмм (переиспользуем память)
            for (int t = 0; t < numThreads; ++t) {
                std::fill(histThreads[t].begin(), histThreads[t].end(), 0);
            }
            std::fill(globalHist.begin(), globalHist.end(), 0);
        }

        // Оптимизированное копирование для тривиально копируемых типов
        template<typename T = ValueType>
        void copyValues(T* dest, const T* src, size_t n, enable_if_trivially_copyable<T>* = nullptr) {
            std::memcpy(dest, src, n * sizeof(T));
        }

        // Копирование для нетривиально копируемых типов
        template<typename T = ValueType>
        void copyValues(T* dest, const T* src, size_t n, enable_if_not_trivially_copyable<T>* = nullptr) {
            for (size_t i = 0; i < n; ++i) {
                dest[i] = src[i];
            }
        }

    public:
        OptimizedRadixSorter() : numThreads(0), currentCapacity(0), initialized(false) {
#pragma omp parallel
            {
#pragma omp single
                numThreads = omp_get_num_threads();
            }
            if (numThreads == 0) numThreads = 1;

            initializeThreadData(numThreads);
        }

        // Конструктор с предварительным выделением памяти
        explicit OptimizedRadixSorter(size_t maxSize) : numThreads(0), currentCapacity(0), initialized(false) {
#pragma omp parallel
            {
#pragma omp single
                numThreads = omp_get_num_threads();
            }
            if (numThreads == 0) numThreads = 1;

            initializeThreadData(numThreads);
            ensureCapacity(maxSize);
        }

        // Основной метод сортировки
        void sort(int* keys, ValueType* values, size_t n) {
            if (n < 2) return;

            ensureCapacity(n);
            resetHistograms();

            for (int pass = 0; pass < PASSES; ++pass) {
                int shift = pass * 8;

                // 1. Параллельный подсчёт гистограмм
#pragma omp parallel for schedule(static)
                for (int i = 0; i < n; ++i) {
                    int tid = omp_get_thread_num();
                    int bucket = (keys[i] >> shift) & MASK;
                    histThreads[tid][bucket]++;
                }

                // 2. Объединение гистограмм
                for (int t = 0; t < numThreads; ++t) {
                    const auto& hist = histThreads[t];
                    auto& global = globalHist;
                    for (int b = 0; b < RADIX; ++b) {
                        global[b] += hist[b];
                    }
                }

                // 3. Вычисление префиксных сумм
                std::vector<int> prefix(RADIX);
                int sum = 0;
                for (int b = 0; b < RADIX; ++b) {
                    prefix[b] = sum;
                    sum += globalHist[b];
                }

                // 4. Расчёт стартовых позиций для каждого потока
                std::copy(prefix.begin(), prefix.end(), localStarts[0].begin());

                for (int t = 1; t < numThreads; ++t) {
                    auto& start = localStarts[t];
                    const auto& prevStart = localStarts[t - 1];
                    const auto& prevHist = histThreads[t - 1];

                    for (int b = 0; b < RADIX; ++b) {
                        start[b] = prevStart[b] + prevHist[b];
                    }
                }

                // 5. Копируем стартовые позиции в счётчики
                for (int t = 0; t < numThreads; ++t) {
                    std::copy(localStarts[t].begin(), localStarts[t].end(),
                        localCounters[t].begin());
                }

                // 6. Параллельная разноска элементов
#pragma omp parallel for schedule(static)
                for (int i = 0; i < n; ++i) {
                    int tid = omp_get_thread_num();
                    int bucket = (keys[i] >> shift) & MASK;
                    int pos = localCounters[tid][bucket]++;

                    keyBuffer[pos] = keys[i];
                    valueBuffer[pos] = values[i];
                }

                // 7. Копирование обратно в исходные массивы
#pragma omp parallel for schedule(static)
                for (int i = 0; i < n; ++i) {
                    keys[i] = keyBuffer[i];
                    values[i] = valueBuffer[i];
                }

                // 8. Сброс гистограмм
                for (int t = 0; t < numThreads; ++t) {
                    std::fill(histThreads[t].begin(), histThreads[t].end(), 0);
                }
                std::fill(globalHist.begin(), globalHist.end(), 0);
            }
        }

        // Перегруженный метод для работы с std::vector
        void sort(std::vector<int>& keys, std::vector<ValueType>& values) {
            if (keys.size() != values.size()) {
                throw std::runtime_error("Key and value arrays must have same size");
            }
            sort(keys.data(), values.data(), keys.size());
        }

        // Метод для сортировки с указанием диапазона
        void sort(int* keys, ValueType* values, size_t start, size_t end) {
            if (start >= end) return;
            sort(keys + start, values + start, end - start);
        }

        void releaseMemory() {
            std::vector<int>().swap(keyBuffer);
            std::vector<ValueType>().swap(valueBuffer);
            std::vector<int>().swap(globalHist);
            std::vector<std::vector<int>>().swap(histThreads);
            std::vector<std::vector<int>>().swap(localStarts);
            std::vector<std::vector<int>>().swap(localCounters);
            currentCapacity = 0;
            initialized = false;
        }

        size_t getMemoryUsage() const {
            size_t total = keyBuffer.capacity() * sizeof(int);
            total += valueBuffer.capacity() * sizeof(ValueType);
            total += globalHist.capacity() * sizeof(int);
            for (const auto& hist : histThreads) {
                total += hist.capacity() * sizeof(int);
            }
            for (const auto& start : localStarts) {
                total += start.capacity() * sizeof(int);
            }
            for (const auto& counter : localCounters) {
                total += counter.capacity() * sizeof(int);
            }
            return total;
        }

        // Получить число потоков
        int getNumThreads() const { return numThreads; }

        // Установить число потоков (должно быть вызвано ДО первой сортировки)
        void setNumThreads(int threads) {
            if (!initialized) {
                numThreads = threads;
                initializeThreadData(numThreads);
            }
        }
    };

    // Частичная специализация для указателей (если нужно хранить указатели)
    template<typename ValueType>
    class OptimizedRadixSorter<ValueType*> {
    private:
        std::vector<int> keyBuffer;
        std::vector<ValueType*> valueBuffer;
        // ... остальная реализация аналогична
    };

    // Примеры использования
    struct ComplexData {
        double x;
        double y;
        int id;

        ComplexData() : x(0), y(0), id(0) {}
        ComplexData(double _x, double _y, int _id) : x(_x), y(_y), id(_id) {}

        bool operator==(const ComplexData& other) const {
            return x == other.x && y == other.y && id == other.id;
        }
    };

}//VM2D

#endif

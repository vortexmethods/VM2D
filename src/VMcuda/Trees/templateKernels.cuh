/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: temploateKernels.cuh                                             |
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
\brief Заголовки и реализации шаблонных функций работы с деревьями
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef TEMPLATEKERNELS_H    
#define TEMPLATEKERNELS_H
	
#include "operations.cuh"


namespace BHcu
{
    template <int order>
    __global__
        __launch_bounds__(THREADSforces, FACTORforces)
        void treeVorticesToPointsCalculationKernel(
            const int nnodesd,   	//количество узлов дерева вихрей (листья-вихри + внутренние узлы, и это число расширено до общего числа тредов)
            const int nbodiesd,   	//количество вихрей
            const double itolsqd, 	//theta^2 в критерии близости
            const double epssqd,  	//радиус вихря (Рэнкина, надо доделать обработку вихрей Ламба)
            const int2* __restrict Mchildd,  //массив индексов потомков для внутренних ячеек дерева вихрей
            const double2* __restrict momsd, //массив мультипольных моментов внутренних ячеек дерева вихрей
            const double3* __restrict vtxd,  //массив вихрей (r.x, r.y, G), в произвольном несортированном порядке
            const int* __restrict MmortonCodesIdxd, //перенумерация исходных индексов для отсортированного массива вихрей, т.е. origVortexIndex = MmortonCodesIdxd[mortonVortexIndex]
            const double2* __restrict Mposd,   //массив центров ячеек дерева 
            const int* __restrict MindexSortd, //перенумерация внутренних ячеек дерева, karrasOrder = MindexSortd[reversed@burtscherOrder], по возрастанию уровня level
            const int* __restrict MindexSortTd,//перенумерация внутренних ячеек дерева, reversed@burtscherOrder = MindexSortTd[karrasOrder]
            const double4* __restrict Mlowerupperd,       //массив координат левых нижних и правых верхних углов внутренних ячеек дерева вихрей    

            const int npointsd,     			//количество точек наблюдения
            const double3* __restrict pointsd,  //массив точек наблюдения в формате (x,y,G), параметр G не используется
            const int* __restrict MmortonCodesIdxPointsd, //перенумерация исходных индексов для отсортированного массива точек наблюдения, т.е. origPointIndex = MmortonCodesIdxPointsd[mortonPointIndex]
            double2* __restrict veld,					  //массив для результатов вычисления скоростей

            bool calcEpsAst,                    //признак вычисления eps*
            double* __restrict epsast          //массив eps*
        )
    {
        register int base = threadIdx.x / WARPSIZE; //номер варпа, в который входит данный тред
        register int sbase = base * WARPSIZE;       //номер первой (точнее, "нулевой") нити в том варпе, в который входит данный тред (ее индекс делится нацело на 32)    	
        register int j = base * MAXDEPTH;           //смещения в массивах pos и node для всех нитей данного варпа
        register int diff = threadIdx.x - sbase;    //номер нити в своем варпе (0 <= diff <= 31)
        register int depth;                         //позиция верхушки стека для всех тредов данного варпа (стек пуст, если depth >= j)

        int npointsUp = ((npointsd + WARPSIZE - 1) / WARPSIZE) * WARPSIZE; //округляем число точек наблюдения до целого числа варпов (до кратного 32)

        register int k;       //индекс обрабатываемой точки наблюдения в порядке mortonOrder
        register int nd;      // индекс родительской ячейки дерева, которую обходим, и которая находится на верхушке стека; 
        register int2 chBoth; //индексы обоих потомков узла nd
        register int pd;      // 0 или 1 --- какого потомка ячейки nd обходим (левого или правого)   

        register const double2* mom; // указатель на набор мультипольных моментов влияющей ячейки в глобальной памяти

        // для вычисления квадратов расстояний до трех ближайших вихрей
        register double d_1 = 1e+5, d_2 = 1e+5, d_3 = 1e+5, dst23, dst12;

        register double3 pointTemp;//временная регистровая переменная для чтения вихря (x,y,G)

        register int indexOfPoint; //истинный индекс точки наблюдения	
        register int srtT;         //истинный индекс внутреннего узла дерева

        register double2 p;        //координаты точки наблюдения

        register double2 v;        //результат расчета скорости в точке наблюдения

        bool isVortex;             //признак того, что обрабатываемая вершина - лист
        register double2 ps;       //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
        register double gm;        //циркуляция вихря, если обрабатываемая вершина --- лист, или 0-й мультипольный момент (суммарная циркуляция вихрей) если обрабатывается внутренняя ячейка дерева
        register double sumSide2;  //габариты обрабатываемой ячейки
        register double2 dr;       //радиус-вектор из центра влияющей ячейки в точку наблюдения    
        register double r2;        //квадрат модуля предыдущего

        __shared__ volatile int pos[MAXDEPTH * THREADSforces / WARPSIZE], node[MAXDEPTH * THREADSforces / WARPSIZE];
        __shared__ double2 momSh[(THREADSforces / WARPSIZE) * orderAlignment]; //массив для подгрузки мультипольных моментов влияющей ячейки 

        __syncthreads();       //синхронизация всех нитей в блоке
        __threadfence_block(); //нулевая нить ждет готовности всех нитей своего блока (по идее, здесь это ни к чему!)

        // iterate over all bodies assigned to thread

        // число блоков точно равно аппаратному числу мультипроцессоров
        // размер блока кратен размеру варпа
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsUp; k += blockDim.x * gridDim.x)
            //0-й блок - точки наблюдения          0 ...   blockDim - 1
            //1-й блок - точки наблюдения   blockDim ... 2*blockDim - 1
            //2-й блок - точки наблюдения 2*blockDim ... 3*blockDim - 1
            //...
            //последний блок - точки                 ... npointsd - 1 + еще несколько "пустых тредов"	
        {

            //только для реальных тел (без "пустых тредов") 
            if (k < npointsd)
            {
                d_1 = d_2 = d_3 = 1e+5;

                indexOfPoint = MmortonCodesIdxPointsd[k];   //истинный индекс точки наблюдения
                pointTemp = pointsd[indexOfPoint];          //координаты точки наблюдения -> во временную регистровую переменную
                p = make_double2(pointTemp.x, pointTemp.y); //они же, но в формате double2 для удобства дальнейшей арифметики

                v.x = v.y = 0; //обнуление результата
            }//if (k < npointsd) 
            //конец кода для "непустых тредов" 

            // initialize iteration stack, i.e., push root node onto stack

            depth = j; //позиция в стеке, с которой надо начать работу

            if (sbase == threadIdx.x) //нулевая нить варпа (можно написать diff==0)
            {
                pos[j] = 0;           //значит, начинаем обрабатывать с левого потомка
                node[j] = nnodesd - 1;//корень дерева в порядке burtscherOrder (последний элемент массива)
            }// if sbase ==

            do
            {
                // stack is not empty
                pd = pos[depth];       //если pd==0, то обрабатываем левого потомка, если pd==1 --- правого потомка
                nd = node[depth];	   //узел, потомков которого обрабаываем

                //karrasOrder:
                // 0 1 2 ... (nb-2) x (nb+0) (nb+1) (nb+2) ... (nb+(nb-1))
                // ----------------   -----------------------------------
                //      cells                         bodies

                //burtscherOrder
                // 0 1 2 ... (nb-1) x x x x (nn-(nb-1)) ... (nn-2) (nn-1)
                // ----------------          ----------------------------
                //      bodies                 sorted and reversed cells

                //в порядке burtscherOrder внутренние узлы пронумерованы подряд с конца: 
                //nd - номер внутреннего узла дерева в порядке burtscherOrder (корень --- (nn-1)-й)
                //(nnodesd - 1) - nd --- номер внутреннего узла в "развернутом порядке" (корень --- 0-й)
                //MindexSortd[(nnodesd - 1) - nd] --- номер внутреннего узла в исходном порядке karrasOrder

                chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]]; //читаем индексы обоих потомков узла nd (в порядке karrasOrder)

                while (pd < 2) //обходим обоих потомков
                {
                    // node on top of stack has more children to process

                    //////////////////////////////////////////////////////////////////////////////////////
                    // Код выполняется всеми тредами варпа, включая "пустые треды"                      //
                    // -------------------------------------------------------------------------------- //
                    // Значения следующих переменных --- одинаковые у всего варпа!!!                    //
                    // chBoth, pd, chd, isVortex, (n, ps, gm, sumSide2) или (n, ps, mom, gm, sumSide2)  //
                    //////////////////////////////////////////////////////////////////////////////////////

                    register int chd = pd * chBoth.y + (1 - pd) * chBoth.x; //индекс обрабатываемого потомка в порядке karrasOrder (если pd == 0, то chBoth.x; если pd == 1, то chBoth.y)
                    ++pd;//для следующего захода перемещаем счетчик pd с 0 на 1

                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    // далее всегда обрабатываем ячейку chd --- это номер внутреннего узла в порядке karrasOrder //
                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    isVortex = (chd >= nbodiesd); //признак того, что обрабатывамая ячейка --- лист (в правой половине дерева в порядке karrasOrder)

                    register int n; //для хранения индекса вихря или индекса ячейки

                    if (isVortex)//если лист
                    {
                        n = chd - nbodiesd;                          //номер вихря в мортоновском порядке; MmortonCodesIdxd[n] --- истинный номер вихря
                        pointTemp = vtxd[MmortonCodesIdxd[n]];       //координаты влияющего вихря -> во временную регистровую переменную 

                        ps = make_double2(pointTemp.x, pointTemp.y); //координаты влияющего вихря
                        gm = pointTemp.z;                            //циркуляция влияющего вихря
                        sumSide2 = 0.0;                     	     //листовая ячейка (вихрь) размера не имеет, т.к. является точкой
                    }//if (isVortex)
                    else//если не лист, а ячейка
                    {
                        srtT = MindexSortTd[chd];   //номер внутреннего узла в "развернутом порядке burtscherOrder" (когда корень --- 0-й), отвечающий узлу chd
                        n = (nnodesd - 1) - srtT;   //номер внутреннего узла в порядке burtscherOrder (когда корень --- последний)
                        ps = Mposd[chd];			//координаты центра внутреннего узла --- влияющей ячейки

                        double4 gab = Mlowerupperd[chd];

                        sumSide2 = (gab.z - gab.x + gab.w - gab.y);
                        sumSide2 *= sumSide2;                    //квадрат суммы габаритов влияющей ячейки
                    }//else if (isVortex)                

                    if (k >= npointsd)
                        sumSide2 = 0.0;

                    dr = p - ps; //радиус-вектор из центра влияющей ячейки в точку наблюдения                
                    r2 = (k < npointsd) ? (dr.x * dr.x + dr.y * dr.y) : 1000000.0; //для непустых тредов --- квадрат расстояния до центра влияющей ячейки

                    // проверка того, что либо влияющий объект - вихрь, либо ячейка, но расположенная достаточно далеко от всех точек наблюдения, обрабатываемых одним варпом                
                    if (isVortex || __all_sync(0xffffffff, (sumSide2 + epssqd) * itolsqd < r2))
                    {
                        if (!isVortex)
                        {
                            mom = momsd + (srtT * orderAlignment);       //указатель на набор мультипольных моментов влияющей ячейки

                            //те треды в варпе, номера которых соответствуют учитываемым моментам, загружают эти моменты влияющей ячейки в shared-memory (для своего варпа внутри блока)
                            if (diff < order)
                                momSh[base * orderAlignment + diff] = mom[diff];

                            gm = momSh[base * orderAlignment + 0].x; //mom[0].x   //нулевой момент (суммарная циркуляция вихрей в ячейке)
                        }//if (!isVortex)

                        //расчет выполняют только "непустые треды"
                        if (k < npointsd)
                        {
                            //если надо найти расстояния до трех ближайших соседей --- пока не важно до отдельных вихрей или до кластеров
                            //todo НАВЕРНОЕ, если влияние от кластера --- то можно не считать, т.к. по определению слишком далеко!
                            if (calcEpsAst)
                            {
                                if ((r2 < d_3) && (r2 > 0))     //если текущий третий сосед должен быть вытолкнут
                                {
                                    dst23 = ::fmin(r2, d_2);    //меньшее из двух чисел --- вновь найденное расстояние и старое второе

                                    d_3 = ::fmax(r2, d_2);      //на третье место --- или вновь вставляемый, или старое второе
                                    //если вновь вставляемое, то dst23 = d_2 (старое второе),
                                    //если старое второе, то dst23 = r2 (вновь найденное)

                                    dst12 = ::fmin(dst23, d_1); //меньшее из чисел --- dst23 и старое первое
                                    d_2 = ::fmax(dst23, d_1);   //на второе место --- или то, что в переменной dst_23, или старое первое

                                    d_1 = dst12; //самое маленькое --- на первую позицию
                                }// if((r2 <...
                            }//if (calcEpsAst)

                            //todo ДОДЕЛАТЬ обработку вихря Ламба
                            register double f = gm / ::fmax(r2, epssqd); //влияние точечного вихря или нулевого момента кластера						
                            v += f * dr;


                            if ((!isVortex) && (order > 1)) //если ячейка, то рассчитываем влияние высших моментов
                            {
                                register double2 cftr; //для расчета коэффициентов theta
                                register double2 th;   //текущее значение theta

                                th = cftr = (r2 ? (1.0 / r2) : 0.0) * dr; //инициализируем их обоих одинаково                                                        
#pragma unroll
                                for (int s = 1; s < order; ++s)
                                {
                                    th = multz(th, cftr);
                                    v += multzA(th, momSh[base * orderAlignment + s]);
                                }//for s

                            }//if ((!isVortex...)

                        }//if (k < npointsd)
                        //конец кода для "непустых тредов"
                    }
                    else //если ячейка, и расположена далеко хотя бы от одной из точек наблюдения в варпе
                    {
                        if (depth < MAXDEPTH * THREADSforces / WARPSIZE) //проверяем, что не вылетели за размер стека
                        {
                            if (pd == 1) //если это была обработка левого потомка, и по нему пришлось идти вниз по дереву
                            {
                                if (sbase == threadIdx.x) //выполняет только нулевой тред в каждом варпе (можно написать diff==0)
                                {
                                    pos[depth] = 1; /*pd*/; //закидываем на стек правого потомка
                                    node[depth] = nd;       //и его родителя
                                    //чтобы потом при развороте стека его рассмотреть
                                }//if sbase

                                ++depth; //сдвигаем верхушку стека вперед, т.к. там появилась новая запись
                            }

                            nd = n; //в качестве базовой ячейки берем ту, для которой принято решение спускаться вниз по дереву
                            pd = 0; //начинаем обход с ее левого потомка

                            chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]]; //подгружаем индексы пары потомков
                        }//if depth <
                    }
                }//while (pd < 2)

                --depth; //когда обошли обоих потомков --- сдвигаем верхушку стека назад, т.к. задание выполнено
            } while (depth >= j); //повторяем процедуру пока стек непустой


            //Сохраняем результат
            if (k < npointsd)
            {
                double2 result = make_double2(-idpid * v.y, idpid * v.x); //реально считали влияние источника, попутно векторно умножаем его на вектор k
                veld[indexOfPoint] = result;

                if (calcEpsAst)
                    epsast[indexOfPoint] = sqrt((d_1 + d_2 + d_3) / 3);
            }//if (k < ndointsd)
        }//for k
    }//treeVorticesToPointsCalculationKernel(...)





    template <int order>
    __global__
        __launch_bounds__(WARPSIZE, FACTORrhs)
        void treeRhsCalculationKernel
        (
            const int nnodesd,      //количество узлов дерева вихрей (листья-вихри + внутренние узлы, и это число расширено до общего числа тредов)
            const int nbodiesd,     //количество вихрей
            const double itolsqd,   //theta^2 в критерии близости
            const int2* __restrict Mchildd,  //массив индексов потомков для внутренних ячеек дерева вихрей
            const double2* __restrict momsd, //массив мультипольных моментов внутренних ячеек дерева вихрей
            const double3* __restrict vtxd,  //массив вихрей (r.x, r.y, G), в произвольном несортированном порядке
            const int* __restrict MmortonCodesIdxd, //перенумерация исходных индексов для отсортированного массива вихрей, т.е. origVortexIndex = MmortonCodesIdxd[mortonVortexIndex]

            const double2* __restrict Mposd,   //массив центров ячеек дерева 
            const int* __restrict MindexSortd, //перенумерация внутренних ячеек дерева, karrasOrder = MindexSortd[reversed@burtscherOrder], по возрастанию уровня level
            const int* __restrict MindexSortTd,//перенумерация внутренних ячеек дерева, reversed@burtscherOrder = MindexSortTd[karrasOrder]
            const double4* __restrict Mlowerupperd, //массив координат левых нижних и правых верхних углов внутренних ячеек дерева вихрей

            const int npointsd,                  //количество панелей на теле
            const double4* __restrict dev_ptr_pt, //массив панелей на теле в формате (beg.x, beg.y, end.x, end.y)
            const double2* __restrict pointsd,   //массив центров панелей в формате (xc, yc)

            const int* __restrict MmortonCodesIdxPointsd,//перенумерация исходных индексов для отсортированного массива панелей, т.е. origPointIndex = MmortonCodesIdxPointsd[mortonPointIndex]

            double* __restrict veld,               //массив для результатов вычисления правых частей (константная проекционная функция)
            double* __restrict vellind            //массив для результатов вычисления правых частей (линейная проекционная функция)

        )
    {
        const int WIDTH5 = min(blockDim.x, WARPSIZE); // варп или усеченный варп, если число тредов в блоке меньше 32

        register int base = threadIdx.x / WIDTH5; //номер варпа или усеченного варпа, в который входит данный тред
        register int sbase = base * WIDTH5;       //номер первой (точнее, "нулевой") нити в том (усеченном) варпе, в который входит данный тред (ее индекс делится на WIDTH5)    	
        register int j = base * MAXDEPTH;         //смещения в массивах pos и node для всех нитей данного (усеченного) варпа
        register int diff = threadIdx.x - sbase;  //номер нити в своем (усеченном) варпе (0 <= diff <= WIDTH5-1)
        register int depth;                       //позиция верхушки стека для всех тредов данного (усеченного) варпа (стек пуст, если depth >= j)

        int npointsUp = ((npointsd + WIDTH5 - 1) / WIDTH5) * WIDTH5;  //округляем число точек наблюдения до целого числа (усеченных) варпов (до кратного WIDTH5)

        register int k;       //индекс обрабатываемой панели наблюдения в порядке mortonOrder
        register int nd;      //индекс родительской ячейки дерева, которую обходим, и которая находится на верхушке стека; 
        register int2 chBoth; //индексы обоих потомков узла nd
        register int pd;      //0 или 1 --- какого потомка ячейки nd обходим (левого или правого)

        register const double2* mom; //указатель на набор мультипольных моментов влияющей ячейки в глобальной памяти

        register double3 pointTemp;  //временная регистровая переменная для чтения вихря (x,y,G)
        register int indexOfPoint;   //истинный индекс панели наблюдения	
        register int srtT;           //истинный индекс внутреннего узла дерева

        register double2 p;          //координаты центра панели наблюдения
        register double val, vallin; //результат расчета правой части в точке наблюдения

        bool isVortex;               //признак того, что обрабатываемая вершина - лист
        register double2 ps;         //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
        register double gm;          //циркуляция вихря, если обрабатываемая вершина --- лист, для внутренних ячеека дерева не используется

        register double sumSide2;    //габариты обрабатываемой ячейки

        register double2 dr;         //радиус-вектор из центра влияющей ячейки в центр панели наблюдения    
        register double r2;          //квадрат модуля предыдущего

        register double2 beg, end, tau; //характеристики обрабатываемой панели
        register double dlen2, idlen;

        register double2 Eloc[order]; //Коэффициенты локального разложения --- храним на регистрах

        __shared__ volatile int pos[MAXDEPTH * maxTHREADSrhs / WARPSIZE], node[MAXDEPTH * maxTHREADSrhs / WARPSIZE];
        __shared__ double2 momSh[std::max((maxTHREADSrhs / WARPSIZE), 1) * orderAlignment];

        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        // число блоков точно равно аппаратному числу мультипроцессоров
        // размер блока кратен размеру варпа или менее одного варпа (равен WIDTH5)
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsUp; k += blockDim.x * gridDim.x)
            //0-й блок - панели наблюдения          0 ...   blockDim - 1
            //1-й блок - панели наблюдения   blockDim ... 2*blockDim - 1
            //2-й блок - панели наблюдения 2*blockDim ... 3*blockDim - 1
            //...
            //последний блок - панели                 ... npointsd - 1 + еще несколько "пустых тредов"	
        {

            //только для реальных панелей (без "пустых тредов") 
            if (k < npointsd)
            {
                indexOfPoint = MmortonCodesIdxPointsd[k]; //истинный индекс панели наблюдения
                p = pointsd[indexOfPoint];                //координаты центра панели наблюдения

                double4 pnl = dev_ptr_pt[indexOfPoint]; //начало и конец панели наблюдения
                beg = make_double2(pnl.x, pnl.y);   //отдельно начало
                end = make_double2(pnl.z, pnl.w);   //отдельно конец

                dlen2 = (end.x - beg.x) * (end.x - beg.x) + (end.y - beg.y) * (end.y - beg.y); //квадрат длины
                idlen = rsqrt(dlen2);         //обратная длина панели
                tau = (end - beg) * idlen;    //вектор касательной к панели, направленный от начала к концу

                val = 0;    //обнуление результата
                vallin = 0; //обнуление результата
            }//if (k < npointsd) 
            //конец кода для "непустых тредов"

            // initialize iteration stack, i.e., push root node onto stack

            depth = j; //позиция в стеке, с которой надо начать работу

            if (sbase == threadIdx.x)  //нулевая нить варпа (можно написать diff==0)
            {
                pos[j] = 0;            //значит, начинаем обрабатывать с левого потомка
                node[j] = nnodesd - 1; //корень дерева в порядке burtscherOrder (последний элемент массива)
            }// if sbase ==

            do
            {
                // stack is not empty
                pd = pos[depth];       //если pd==0, то обрабатываем левого потомка, если pd==1 --- правого потомка
                nd = node[depth];	   //узел, потомков которого обрабаываем

                //karrasOrder:
                // 0 1 2 ... (nb-2) x (nb+0) (nb+1) (nb+2) ... (nb+(nb-1))
                // ----------------   -----------------------------------
                //      cells                         bodies

                //burtscherOrder
                // 0 1 2 ... (nb-1) x x x x (nn-(nb-1)) ... (nn-2) (nn-1)
                // ----------------          ----------------------------
                //      bodies                 sorted and reversed cells

                //в порядке burtscherOrder внутренние узлы пронумерованы подряд с конца: 
                //nd - номер внутреннего узла дерева в порядке burtscherOrder (корень --- (nn-1)-й)
                //(nnodesd - 1) - nd --- номер внутреннего узла в "развернутом порядке" (корень --- 0-й)
                //MindexSortd[(nnodesd - 1) - nd] --- номер внутреннего узла в исходном порядке karrasOrder

                chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]]; //читаем индексы обоих потомков узла nd (в порядке karrasOrder)

                while (pd < 2)
                {
                    // node on top of stack has more children to process

                    //////////////////////////////////////////////////////////////////////////////////////
                    // Код выполняется всеми тредами варпа, включая "пустые треды"                      //
                    // -------------------------------------------------------------------------------- //
                    // Значения следующих переменных --- одинаковые у всего варпа!!!                    //
                    // chBoth, pd, chd, isVortex, (n, ps, gm, sumSide2) или (n, ps, mom, gm, sumSide2)  //
                    //////////////////////////////////////////////////////////////////////////////////////

                    register int chd = pd * chBoth.y + (1 - pd) * chBoth.x; //индекс обрабатываемого потомка в порядке karrasOrder (если pd == 0, то chBoth.x; если pd == 1, то chBoth.y)
                    ++pd;//для следующего захода перемещаем счетчик pd с 0 на 1

                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    // далее всегда обрабатываем ячейку chd --- это номер внутреннего узла в порядке karrasOrder //
                    ///////////////////////////////////////////////////////////////////////////////////////////////
                    isVortex = (chd >= nbodiesd); //признак того, что обрабатывамая ячейка --- лист (в правой половине дерева в порядке karrasOrder)

                    register int n; //для хранения индекса вихря или индекса ячейки

                    if (isVortex)//если лист
                    {
                        n = chd - nbodiesd;                         //номер вихря в мортоновском порядке; MmortonCodesIdxd[n] --- истинный номер вихря
                        pointTemp = vtxd[MmortonCodesIdxd[n]];      //координаты влияющего вихря -> во временную регистровую переменную
                        ps = make_double2(pointTemp.x, pointTemp.y);//координаты влияющего вихря
                        gm = pointTemp.z;                           //циркуляция влияющего вихря
                        sumSide2 = 0;                               //листовая ячейка (вихрь) размера не имеет, т.к. является точкой
                    }//if (isVortex)
                    else//если не лист, а ячейка
                    {
                        srtT = MindexSortTd[chd];              //номер внутреннего узла в "развернутом порядке burtscherOrder" (когда корень --- 0-й), отвечающий узлу chd
                        n = (nnodesd - 1) - srtT;              //номер внутреннего узла в порядке burtscherOrder (когда корень --- последний)
                        ps = Mposd[chd];                       //координаты центра внутреннего узла --- влияющей ячейки

                        double4 gab = Mlowerupperd[chd];          //габарит влияющей ячейки
                        sumSide2 = gab.z - gab.x + gab.w - gab.y; //сумма габаритов (полупериметр)
                        sumSide2 *= sumSide2;                     //квадрат суммы габаритов влияющей ячейки              
                    }
                    if (k >= npointsd)
                        sumSide2 = dlen2 = 0.0;

                    dr = p - ps;        //радиус-вектор из центра влияющей ячейки в центр панели наблюдения   
                    r2 = (k < npointsd) ? (dr.x * dr.x + dr.y * dr.y) : 1000000.0;  //для непустых тредов --- квадрат расстояния до центра влияющей ячейки

                    // проверка того, что либо влияющий объект - вихрь, либо ячейка, но расположенная достаточно далеко от всех панелей наблюдения, обрабатываемых одним варпом   
                    if (isVortex || __all_sync(0xffffffff, (sumSide2 + dlen2) * itolsqd < r2))
                    {
                        if (!isVortex)
                        {
                            mom = momsd + (srtT * orderAlignment);                  //указатель на набор мультипольных моментов влияющей ячейки

                            //те треды в варпе, номера которых соответствуют учитываемым моментам, загружают эти моменты влияющей ячейки в shared-memory (для своего варпа внутри блока)
                            if (diff < order)
                                momSh[base * (orderAlignment)+diff] = mom[diff];
                            //__threadfence_block(); //синхронизация внутри блока; не нужна, т.к. варп и так работает синхронно, а требуется именно синхронность варпа!

                        }//if (!isVortex)              

                        //расчет выполняют только "непустые треды"
                        if (k < npointsd)
                        {
                            //влияние точечного вихря на панели
                            if (isVortex)
                            {
                                double2 ss = ps - beg;
                                double2 pp = ps - end;

                                double alpha = atan2(pp.x * ss.y - pp.y * ss.x, pp.x * ss.x + pp.y * ss.y);
                                double tempVel = gm * alpha;

                                if (r2 > 1e-20)
                                    val -= tempVel;
                                //else
                                //    val -= 0.5 /idpi * gm;                       

                                //расчет с линейной проекционной функцией, если используем схему T1
                                if (vellind != nullptr)
                                {
                                    double txx = tau.x * tau.x;
                                    double txy = tau.x * tau.y;
                                    double tyy = tau.y * tau.y;

                                    double2 u1;
                                    u1.x = (pp.x + ss.x) * (txx - tyy) + 2.0 * (pp.y + ss.y) * txy;
                                    u1.y = (pp.y + ss.y) * (tyy - txx) + 2.0 * (pp.x + ss.x) * txy;

                                    double lambda = 0.5 * log((ss.x * ss.x + ss.y * ss.y) / (pp.x * pp.x + pp.y * pp.y));

                                    double tempVelLin = gm * (alpha * (u1.x * tau.x + u1.y * tau.y) + lambda * (-u1.y * tau.x + u1.x * tau.y));

                                    vallin -= 0.5 * idlen * tempVelLin;
                                }//if (vellind != nullptr)
                            }//if (isVortex)

                            //если влияние считаем от ячейки, содержащей несколько вихрей
                            else
                            {

                                //рассчитываем коэффициенты локального разложения
                                double2 theta = (p - ps) / r2;

#pragma unroll
                                for (int q = 0; q < order; ++q)
                                    Eloc[q] = make_double2(0.0, 0.0);

#pragma unroll
                                for (int q = 0; q < order; ++q)
                                {
                                    for (int s = q; s >= 0; --s)
                                        Eloc[s] += (((2 * q - s) & 1) ? -ifac[q - s + 1] : ifac[q - s + 1]) * multzA(theta, momSh[base * (orderAlignment)+q - s]);
                                    theta = ((q + 1) / r2) * multz(theta, p - ps);
                                }

                                double2 v = 0.5 * Eloc[0];
                                double2 vL{ 0.0, 0.0 };

                                double2 rPan = end - beg;

                                double2 kp;
                                double2 mulP = kp = rPan;

                                double2 taudL = (0.5 / dlen2) * rPan;
                                double2 taudLc = taudL;

#pragma unroll
                                for (int k = 1; k < order; ++k)
                                {
                                    mulP = multz(mulP, kp);
                                    taudL *= 0.5;

                                    if (!(k & 1))
                                        v += ifac[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                    else
                                        if (vellind != nullptr)
                                            vL += ifac[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
                                }

                                val += 2.0 * (-v.y * rPan.x + v.x * rPan.y);

                                if (vellind != nullptr)
                                    vallin += 2.0 * (-vL.y * rPan.x + vL.x * rPan.y);
                            }//рассчитали влияние от ячейки целиком по мультипольной схеме

                        }//if (k < npointsd)
                        //конец кода для "непустых тредов"
                    }
                    else //если ячейка, и расположена далеко хотя бы от одной из точек наблюдения в варпе
                    {
                        if (depth < MAXDEPTH * maxTHREADSrhs / WARPSIZE) //проверяем, что не вылетели за размер стека
                        {
                            if (pd == 1) //если это была обработка левого потомка, и по нему пришлось идти вниз по дереву
                            {
                                if (sbase == threadIdx.x) //выполняет только нулевой тред в каждом варпе (можно написать diff==0)
                                {
                                    pos[depth] = 1; /*pd*/; //закидываем на стек правого потомка
                                    node[depth] = nd;       //и его родителя
                                    //чтобы потом при развороте стека его рассмотреть
                                }

                                ++depth; //сдвигаем верхушку стека вперед, т.к. там появилась новая запись
                            }
                            nd = n; //в качестве базовой ячейки берем ту, для которой принято решение спускаться вниз по дереву
                            pd = 0; //начинаем обход с ее левого потомка

                            chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];  //подгружаем индексы пары потомков
                        }//if depth <
                    }
                }//while (pd < 2)

                --depth; //когда обошли обоих потомков --- сдвигаем верхушку стека назад, т.к. задание выполнено

            } while (depth >= j); //повторяем процедуру пока стек непустой

            //Сохраняем результат
            if (k < npointsd)
            {
                register double cf = idpid * idlen;
                veld[indexOfPoint] = cf * val;

                if (vellind != nullptr)
                    vellind[indexOfPoint] = cf * vallin;
            }//if (k < ndointsd)        
        }//for k
    }//treeRhsCalculationKernel(...)





}





	
#endif
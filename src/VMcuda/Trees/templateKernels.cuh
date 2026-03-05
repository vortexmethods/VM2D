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
    template <int order, int THREADSforces>
    __global__
        __launch_bounds__(THREADSforces, FACTORforces)
        void treeVorticesToPointsCalculationKernel(
            const int nnodesd,   	//количество узлов дерева вихрей (листья-вихри + внутренние узлы, и это число расширено до общего числа тредов)
            const int nbodiesd,   	//количество вихрей
            const double itolsqd, 	//theta^2 в критерии близости
            const int2* __restrict Mchildd,  //массив индексов потомков для внутренних ячеек дерева вихрей
            const double2* __restrict momsd, //массив мультипольных моментов внутренних ячеек дерева вихрей
            const double4* __restrict vtxd,  //массив вихрей (r.x, r.y, G, Sigma), в произвольном несортированном порядке
            const int* __restrict MmortonCodesIdxd, //перенумерация исходных индексов для отсортированного массива вихрей, т.е. origVortexIndex = MmortonCodesIdxd[mortonVortexIndex]
            const double2* __restrict Mposd,   //массив центров ячеек дерева 
            const int* __restrict MindexSortd, //перенумерация внутренних ячеек дерева, karrasOrder = MindexSortd[reversed@burtscherOrder], по возрастанию уровня level
            const int* __restrict MindexSortTd,//перенумерация внутренних ячеек дерева, reversed@burtscherOrder = MindexSortTd[karrasOrder]
            const double4* __restrict Mlowerupperd,       //массив координат левых нижних и правых верхних углов внутренних ячеек дерева вихрей    

            const int npointsd,     			//количество точек наблюдения
            const double4* __restrict pointsd,  //массив точек наблюдения в формате (x,y,G,Sigma), параметры G и Sigma не используются
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

        register double4 pointTemp;//временная регистровая переменная для чтения вихря (x,y,G,S)

        register int indexOfPoint; //истинный индекс точки наблюдения	
        register int srtT;         //истинный индекс внутреннего узла дерева

        register double2 p;        //координаты точки наблюдения
        register double sgm2Contr;
        register double2 v;        //результат расчета скорости в точке наблюдения

        bool isVortex;             //признак того, что обрабатываемая вершина - лист
        register double2 ps;       //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
        register double gm;        //циркуляция вихря, если обрабатываемая вершина --- лист, или 0-й мультипольный момент (суммарная циркуляция вихрей) если обрабатывается внутренняя ячейка дерева
        register double sgm2;
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
                sgm2Contr = sqr(pointTemp.w);

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
                        sgm2 = sqr(pointTemp.w);
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
                    if (isVortex || __all_sync(0xffffffff, (sumSide2 + sgm2Contr) * itolsqd < r2))
                    {
                        if (!isVortex)
                        {
                            mom = momsd + (srtT * orderAlignment);       //указатель на набор мультипольных моментов влияющей ячейки

                            //те треды в варпе, номера которых соответствуют учитываемым моментам, загружают эти моменты влияющей ячейки в shared-memory (для своего варпа внутри блока)
                            if (diff < order)
                                momSh[base * orderAlignment + diff] = mom[diff];

                            gm = momSh[base * orderAlignment + 0].x; //mom[0].x   //нулевой момент (суммарная циркуляция вихрей в ячейке)
                            sgm2 = 0;
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
                            register double f = gm / ::fmax(r2, sgm2); //влияние точечного вихря или нулевого момента кластера						
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
            const double4* __restrict vtxd,  //массив вихрей (r.x, r.y, G, S), в произвольном несортированном порядке
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

        register double4 pointTemp;  //временная регистровая переменная для чтения вихря (x,y,G,S)
        register int indexOfPoint;   //истинный индекс панели наблюдения	
        register int srtT;           //истинный индекс внутреннего узла дерева

        register double2 p;          //координаты центра панели наблюдения
        register double val, vallin; //результат расчета правой части в точке наблюдения

        bool isVortex;               //признак того, что обрабатываемая вершина - лист
        register double2 ps;         //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
        register double gm;          //циркуляция вихря, если обрабатываемая вершина --- лист, для внутренних ячеека дерева не используется
        //register double sgm;
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
                        //sgm = pointTemp.w;                           
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
                                        Eloc[s] += (((2 * q - s) & 1) ? -ifacGPU[q - s + 1] : ifacGPU[q - s + 1]) * multzA(theta, momSh[base * (orderAlignment)+q - s]);
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
                                        v += ifacGPU[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                    else
                                        if (vellind != nullptr)
                                            vL += ifacGPU[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
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

    /******************************************************************************/
    /***************  Auxilary functions for SLAE *********************************/
    /******************************************************************************/

    /// \brief Вспомогательная функция, которая определяет, находится ли панель itI на контуре после панели itJ
    ///
    /// \param[in] itI константная ссылка на обертку для второй ("правой") панели
    /// \param[in] itJ константная ссылка на обертку для первой ("левой") панели
    __device__ __forceinline bool isEqual(const double2& paniBeg, const double2& panjEnd)
    {
        //return ((paniBeg.x - panjEnd.x) * (paniBeg.x - panjEnd.x) + (paniBeg.y - panjEnd.y) * (paniBeg.y - panjEnd.y) < 1e-20);
        return (fabs(paniBeg.x - panjEnd.x) + fabs(paniBeg.y - panjEnd.y) < 1e-10);
    }

    /// \brief Вспомогательная функция вычисления угла между векторами (в диапазоне (-pi...pi]) (со знаком, поворот от первого ко второму)
    ///
    /// \param[in] p константная ссылка на первый вектор
    /// \param[in] s константная ссылка на второй вектор 
    __device__ __forceinline double Alpha(const double2& p, const double2& s)
    {
        return atan2(p.x * s.y - p.y * s.x, p.x * s.x + p.y * s.y);
    }

    /// \brief  Вспомогательная функция вычисления логарифма отношения норм векторов
    /// 
    /// \param[in] p константная ссылка на первый вектор
    /// \param[in] s константная ссылка на второй вектор 
    __device__ __forceinline double Lambda(const double2& p, const double2& s)
    {
        return 0.5 * log((s.x * s.x + s.y * s.y) / (p.x * p.x + p.y * p.y));
    }


    /// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
    __device__ __forceinline double2 Omega(const double2& a, const double2& b, const double2& c)
    {
        double adotb = a.x * b.x + a.y * b.y;
        double acrossb = a.x * b.y - a.y * b.x;
        return double2{ adotb * c.x - acrossb * c.y, adotb * c.y + acrossb * c.x };
    };



    template <int order>
    __global__
        __launch_bounds__(THREADSslae, FACTORslae)
        void treeMatrToVecZeroIterKernel12
        (
            const int nnodesd,  // Число узлов дерева (по вихрям)
            const int nbodiesd, // Число вихрей в пелене

            const double itolsqd, // 1/theta^2
            const int2* __restrict Mchildd, //Массив потомков узлов дерева

            const double2* __restrict momsd, //моменты узлов дерева
            const double* __restrict rpnl,  //начала и концы влияющих панелей
            const double* dev_pnl,

            const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

            const double2* __restrict Mposd, //массив координат центров ячеек дерева
            const int* __restrict MindexSortd, //порядок сортировки внутренних узлов дерева
            const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
            const double4* __restrict Mlowerupperd,//координаты нижнего левого и верхнего правого угла ячейки дерева
            scheme_T schemeType,

            const int npointsd,                  //количество панелей на теле
            const double* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)

            double2* __restrict Ed,                //коэффициенты локального разложения
            const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей

            double* __restrict veld,               //куда сохранять ответ
            double* __restrict vellind,            //куда сохранять ответ

            int* __restrict closeCellsPfl,
            int* __restrict ClosePrefixSuml,
            int* __restrict farCellsPfl,
            int* __restrict FarPrefixSuml,
            double2* __restrict i00save,
            double2* __restrict i01save,
            double2* __restrict i10save,
            double2* __restrict i11save, int iter
        )
    {
        const int WIDTH5 = min(blockDim.x, WARPSIZE);

        register int base = threadIdx.x / WIDTH5;
        register int sbase = base * WIDTH5;
        register int j = base * MAXDEPTH;
        register int diff = threadIdx.x - sbase;
        register int depth;

        int npointsUp = ((npointsd + WIDTH5 - 1) / WIDTH5) * WIDTH5;

        register int k;
        register int nd;
        register int2 chBoth;
        register int pd;

        register const double2* mom;
        register int indexOfPoint;
        register int srtT;

        register double2 p;
        register double val, vallin;

        bool isVortex;
        register double2 ps;
        register double gm;

        register double sumSide2;
        register double2 dr;
        register double r2;

        register double2 beg, end;
        register double dlen2;
        register double2 tau;
        register double2 di;
        register double dilen, idlen;

        register double4 gab;

        register double2 Eloc[order];

        int closecntr, farcntr;

        __shared__ volatile int pos[MAXDEPTH * maxTHREADSslae / WARPSIZE], node[MAXDEPTH * maxTHREADSslae / WARPSIZE];
        __shared__ double2 momSh[std::max((maxTHREADSslae / WARPSIZE), 1) * orderAlignment];

        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsUp; k += blockDim.x * gridDim.x)
        {
            if (k < npointsd)
            {
                indexOfPoint = MmortonCodesIdxPointsd[k];

                closecntr = 0;
                farcntr = 0;

                const double4 pnl = *(double4*)&(dev_ptr_pt[4 * indexOfPoint]);
                beg = make_double2(pnl.x, pnl.y);
                end = make_double2(pnl.z, pnl.w);
                p = 0.5 * (beg + end);

                di = make_double2(end.x - beg.x, end.y - beg.y);
                dlen2 = di.x * di.x + di.y * di.y;

                idlen = rsqrt(dlen2);
                dilen = 1.0 / idlen;
                tau = (end - beg) * idlen;

                val = 0;
                vallin = 0;
            }

            // initialize iteration stack, i.e., push root node onto stack
            depth = j;

            if (sbase == threadIdx.x)
            {
                pos[j] = 0;
                node[j] = nnodesd - 1;
            }

            do
            {
                // stack is not empty
                pd = pos[depth];
                nd = node[depth];

                chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

                register int mcn;

                while (pd < 2)
                {
                    // node on top of stack has more children to process

                    // load child pointer
                    //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                    register int chd = pd * chBoth.y + (1 - pd) * chBoth.x;

                    ++pd;

                    isVortex = (chd >= nbodiesd);

                    register int n;


                    if (isVortex)
                    {
                        n = chd - nbodiesd;

                        mcn = MmortonCodesIdxd[n];

                        gab = *(double4*)&(rpnl[mcn * 4]);
                        const double2 pnlVec = make_double2(gab.z - gab.x, gab.w - gab.y);
                        ps = double2{ 0.5 * (gab.x + gab.z), 0.5 * (gab.y + gab.w) };

                        double lj = sqrt(pnlVec.x * pnlVec.x + pnlVec.y * pnlVec.y);

                        sumSide2 = fabs(pnlVec.x) + fabs(pnlVec.y);
                        sumSide2 *= sumSide2;
                    }
                    else
                    {
                        srtT = MindexSortTd[chd];
                        n = (nnodesd - 1) - srtT;
                        ps = Mposd[chd];

                        gab = Mlowerupperd[chd];
                        sumSide2 = gab.z - gab.x + gab.w - gab.y;
                        sumSide2 *= sumSide2;
                    }

                    if (k >= npointsd)
                        sumSide2 = dlen2 = 0.0;

                    //ps - положение вихря/кластера
                    //p - центр панели
                    dr = p - ps;
                    r2 = (k < npointsd) ? (dr.x * dr.x + dr.y * dr.y) : 1000000.0;   // compute distance squared               

                    // check if all threads agree that cell is far enough away (or is a body)
                    if (isVortex || __all_sync(0xffffffff, (sumSide2 + dlen2) * itolsqd < r2))
                    {
                        if (!isVortex)
                        {
                            mom = momsd + (srtT * orderAlignment);
                            if (diff < order)
                                momSh[base * (orderAlignment)+diff] = mom[diff];
                        }
                        __threadfence_block();
                        __syncthreads();//added

                        if (k < npointsd)
                        {
                            if (isVortex)
                            {
                                closeCellsPfl[ClosePrefixSuml[indexOfPoint] + closecntr] = chd;

                                double tempVelNew;
                                const double2 ptpanBeg = make_double2(gab.x, gab.y);
                                const double2 ptpanEnd = make_double2(gab.z, gab.w);

                                const double2 ptpanVec = make_double2(ptpanEnd.x - ptpanBeg.x, ptpanEnd.y - ptpanBeg.y);

                                gm = dev_pnl[12 * mcn + 6];//dev_ptr_freeVortexSheet[mcn];
                                double gmlin = 0.0;
                                if (schemeType == scheme_T::linScheme)//if (dev_ptr_freeVortexSheetLin != nullptr)
                                    gmlin = dev_pnl[12 * mcn + 9];//dev_ptr_freeVortexSheetLin[mcn];

                                if (r2 > 1e-20)
                                {
                                    double2 i00;
                                    double2 dj;
                                    double djlen;
                                    double ilenj;
                                    double2 tauj;
                                    double2 p1, s1, p2, s2;
                                    double3 alpha, lambda;
                                    bool condBefore, condAfter;
                                    if (iter == 0)
                                    {
                                        dj = make_double2(ptpanEnd.x - ptpanBeg.x, ptpanEnd.y - ptpanBeg.y);
                                        ilenj = rsqrt(dj.x * dj.x + dj.y * dj.y);
                                        djlen = 1 / ilenj;
                                        tauj = make_double2(ilenj * dj.x, ilenj * dj.y);

                                        p1 = make_double2(end.x - ptpanEnd.x, end.y - ptpanEnd.y);
                                        s1 = make_double2(end.x - ptpanBeg.x, end.y - ptpanBeg.y);
                                        p2 = make_double2(beg.x - ptpanEnd.x, beg.y - ptpanEnd.y);
                                        s2 = make_double2(beg.x - ptpanBeg.x, beg.y - ptpanBeg.y);

                                        condBefore = isEqual(beg, ptpanEnd);
                                        condAfter = isEqual(ptpanBeg, end);

                                        alpha = make_double3(\
                                            condAfter ? 0.0 : Alpha(s2, s1), \
                                            Alpha(s2, p1), \
                                            condBefore ? 0.0 : Alpha(p1, p2) \
                                        );

                                        lambda = make_double3(\
                                            condAfter ? 0.0 : Lambda(s2, s1), \
                                            Lambda(s2, p1), \
                                            condBefore ? 0.0 : Lambda(p1, p2) \
                                        );

                                        double2 v00_0 = Omega(s1, tau, tauj);
                                        double2 v00_1 = Omega(di, tau, tauj);

                                        double2 v00_2 = Omega(p2, tau, tauj);

                                        i00 = make_double2(
                                            ilenj * (alpha.x * v00_0.x - alpha.y * v00_1.x + alpha.z * v00_2.x \
                                                - (lambda.x * v00_0.y - lambda.y * v00_1.y + lambda.z * v00_2.y)),
                                            ilenj * (alpha.x * v00_0.y - alpha.y * v00_1.y + alpha.z * v00_2.y \
                                                + (lambda.x * v00_0.x - lambda.y * v00_1.x + lambda.z * v00_2.x))
                                        );

                                        i00save[ClosePrefixSuml[indexOfPoint] + closecntr] = i00;
                                    }
                                    else
                                        i00 = i00save[ClosePrefixSuml[indexOfPoint] + closecntr];

                                    tempVelNew = -gm * (i00.x * tau.x + i00.y * tau.y);
                                    val -= tempVelNew;

                                    if (schemeType == scheme_T::linScheme)
                                    {
                                        double2 i01, i10, i11;

                                        if (iter == 0)
                                        {
                                            double s1len2 = s1.x * s1.x + s1.y * s1.y;

                                            double2 om1 = Omega(s1, tau, tauj);
                                            double2 om2 = Omega(make_double2(s1.x + p2.x, s1.y + p2.y), tauj, tauj);
                                            double sc = (p1.x + s1.x) * tauj.x + (p1.y + s1.y) * tauj.y;

                                            const double hLj = 0.5 * ilenj;
                                            const double pA = hLj * sc;
                                            const double pB = hLj * s1len2;
                                            const double hLij = hLj * dilen;

                                            double2 v01_0 = make_double2(pA * om1.x - pB * tau.x, pA * om1.y - pB * tau.y);
                                            double2 v01_1 = make_double2(hLij * om2.x, hLij * om2.y);

                                            i01 = make_double2(\
                                                ilenj * ((alpha.x + alpha.z) * v01_0.x - (alpha.y + alpha.z) * v01_1.x\
                                                    - (((lambda.x + lambda.z) * v01_0.y - (lambda.y + lambda.z) * v01_1.y) - 0.5 * dilen * tauj.y)),
                                                ilenj * ((alpha.x + alpha.z) * v01_0.y - (alpha.y + alpha.z) * v01_1.y\
                                                    + (((lambda.x + lambda.z) * v01_0.x - (lambda.y + lambda.z) * v01_1.x) - 0.5 * dilen * tauj.x))
                                            );

                                            i01save[ClosePrefixSuml[indexOfPoint] + closecntr] = i01;

                                            double2 om3 = Omega(make_double2(s1.x + p2.x, s1.y + p2.y), tau, tau);

                                            const double hLi = 0.5 * idlen;
                                            const double pC = hLi * ((s1 + s2) & tau);
                                            const double pD = hLi * s1len2;
                                            const double hLji = hLi * djlen;

                                            double2 v10_0 = make_double2(pC * om1.x - pD * tauj.x, pC * om1.y - pD * tauj.y);
                                            double2 v10_1 = make_double2(hLji * om3.x, hLji * om3.y);

                                            i10 = make_double2(
                                                ilenj * (-(alpha.x + alpha.z) * v10_0.x + alpha.z * v10_1.x \
                                                    - ((-(lambda.x + lambda.z) * v10_0.y + lambda.z * v10_1.y) + 0.5 * djlen * tau.y)),
                                                ilenj * (-(alpha.x + alpha.z) * v10_0.y + alpha.z * v10_1.y \
                                                    + ((-(lambda.x + lambda.z) * v10_0.x + lambda.z * v10_1.x) + 0.5 * djlen * tau.x))
                                            );

                                            //11-01-2026
                                            i10save[ClosePrefixSuml[indexOfPoint] + closecntr] = i10;

                                            double2 om4 = Omega(make_double2(s1.x - 3.0 * p2.x, s1.y - 3.0 * p2.y), tau, tauj);
                                            double2 om5 = Omega(di, tauj, tauj);
                                            double2 om6 = Omega(dj, tau, tau);
                                            double sc4 = s1.x * om4.x + s1.y * om4.y;

                                            const double mn = idlen * ilenj / 12.0;
                                            const double qA = mn * 2.0 * sc4 - 0.25;
                                            const double qB = mn * s1len2;

                                            double2 v11_0 = make_double2(
                                                qA * om1.x - qB * (s1.x - 3 * p2.x),
                                                qA * om1.y - qB * (s1.y - 3 * p2.y)
                                            );

                                            double qC = hLij / 6.0;
                                            double qD = hLji / 6.0;

                                            double2 v11_1 = make_double2(qC * om5.x, qC * om5.y);
                                            double2 v11_2 = make_double2(qD * om6.x, qD * om6.y);

                                            i11 = make_double2(
                                                ilenj * ((alpha.x + alpha.z) * v11_0.x - (alpha.y + alpha.z) * v11_1.x - alpha.z * v11_2.x\
                                                    - ((lambda.x + lambda.z) * v11_0.y - (lambda.y + lambda.z) * v11_1.y - lambda.z * v11_2.y \
                                                        + 1.0 / 12.0 * (djlen * tau.y + dilen * tauj.y - 2.0 * om1.y))),
                                                ilenj * ((alpha.x + alpha.z) * v11_0.y - (alpha.y + alpha.z) * v11_1.y - alpha.z * v11_2.y\
                                                    + ((lambda.x + lambda.z) * v11_0.x - (lambda.y + lambda.z) * v11_1.x - lambda.z * v11_2.x \
                                                        + 1.0 / 12.0 * (djlen * tau.x + dilen * tauj.x - 2.0 * om1.x)))
                                            );

                                            i11save[ClosePrefixSuml[indexOfPoint] + closecntr] = i11;
                                        }
                                        else
                                        {
                                            i01 = i01save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                            i10 = i10save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                            i11 = i11save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                        }
                                        double tempVelNewB = -gmlin * (i01.x * tau.x + i01.y * tau.y);
                                        val -= tempVelNewB;

                                        double tempVelNewC = -gm * (i10.x * tau.x + i10.y * tau.y);
                                        double tempVelNewD = -gmlin * (i11.x * tau.x + i11.y * tau.y);

                                        vallin -= (tempVelNewC + tempVelNewD);
                                    }
                                }//if (posI != posJ)

                                ++closecntr;

                            }
                            else
                            {
                                farCellsPfl[FarPrefixSuml[indexOfPoint] + farcntr] = chd;

                                double2 theta = (p - ps) / r2;
#pragma unroll
                                for (int q = 0; q < order; ++q)
                                    Eloc[q] = make_double2(0.0, 0.0);

#pragma unroll
                                for (int q = 0; q < order; ++q)
                                {
                                    for (int s = q; s >= 0; --s)
                                        Eloc[s] += (((2 * q - s) & 1) ? -ifacGPU[q - s + 1] : ifacGPU[q - s + 1]) * multzA(theta, momSh[base * (orderAlignment)+q - s]);
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
                                        v += ifacGPU[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                    else
                                        if (schemeType == scheme_T::linScheme)
                                            vL += ifacGPU[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
                                }

                                v *= 2.0;
                                if (schemeType == scheme_T::linScheme)
                                    vL *= 2.0;

                                val += (-v.y * rPan.x + v.x * rPan.y);

                                if (schemeType == scheme_T::linScheme)
                                    vallin += (-vL.y * rPan.x + vL.x * rPan.y);

                                ++farcntr;
                            }
                        }
                        __threadfence_block();
                        __syncthreads();//added
                    }
                    else //идем глубже по дереву
                    {
                        if (depth < MAXDEPTH * maxTHREADSslae / WARPSIZE)
                        {
                            if (pd == 1)
                            {
                                if (sbase == threadIdx.x)
                                {
                                    pos[depth] = 1;
                                    node[depth] = nd;
                                }
                                ++depth;
                            }
                            pd = 0;
                            nd = n;

                            chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                        }//if depth<
                    }
                }
                --depth;

            } while (depth >= j);

            // update velocity
            if (k < npointsd)
            {
                register double cf = idpid * idlen;
                veld[indexOfPoint] = cf * val;

                if (schemeType == scheme_T::linScheme)
                    vellind[indexOfPoint] = cf * vallin;
            };
        }
    }//treeMatrToVecZeroIterKernel12(...)

    template <int order>
    __global__
        __launch_bounds__(THREADSslae, FACTORslae)
        void treeMatrToVecNonZeroIterKernel12 //new 
        (
            const int nnodesd,  // Число узлов дерева (по вихрям)
            const int nbodiesd, // Число вихрей в пелене

            const double2* __restrict momsd, //моменты узлов дерева
            const double* __restrict rpnl,  //вихри в пелене
            const double* dev_pnl,
            /*const double* dev_ptr_freeVortexSheet,
            const double* dev_ptr_freeVortexSheetLin,*/

            const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

            const double2* __restrict Mposd, //массив координат центров ячеек дерева
            const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
            scheme_T schemeType,

            const int npointsd,                  //количество панелей на теле
            const double* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)

            double2* __restrict Ed,                //коэффициенты локального разложения
            const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей

            double* __restrict veld,               //куда сохранять ответ
            double* __restrict vellind,            //куда сохранять ответ

            int* __restrict nClosePanelsl,
            int* __restrict nFarCellsl,

            const int* __restrict closeCellsPfl,
            const int* __restrict ClosePrefixSuml,
            const int* __restrict farCellsPfl,
            const int* __restrict FarPrefixSuml,

            const double2* __restrict i00save,
            const double2* __restrict i01save,
            const double2* __restrict i10save,
            const double2* __restrict i11save
        )
    {
        const int WIDTH5 = min(blockDim.x, WARPSIZE);
        register int k, n;
        register int base = threadIdx.x / WIDTH5;
        register int sbase = base * WIDTH5;
        register int diff = threadIdx.x - sbase;
        register double2 p, ps;
        register const double2* mom;
        register double val, vallin;

        register int MmortonCodesIdxd_n = 0;

        register double2 Eloc[order];


        bool calcLin = (schemeType == scheme_T::linScheme);


        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread

        int npointsdUp = ((npointsd + 31) / 32) * 32;

        int indexOfPoint;
        double2 beg, end;
        double2 di, tau;
        double dlen2, idlen;
        __shared__ int nClosePnls;
        __shared__ int nFarCells;

        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsdUp; k += blockDim.x * gridDim.x)
        {
            if (k < npointsd)
            {
                indexOfPoint = MmortonCodesIdxPointsd[k];

                double4 begend = *(double4*)(dev_ptr_pt + 4 * indexOfPoint);

                beg = make_double2(begend.x, begend.y);
                end = make_double2(begend.z, begend.w);
                p = 0.5 * (beg + end);

                di = make_double2(end.x - beg.x, end.y - beg.y);
                dlen2 = di.x * di.x + di.y * di.y;

                idlen = rhypot(di.x, di.y);

                tau = (end - beg) * idlen;

                val = 0;
                vallin = 0;

                nClosePnls = nClosePanelsl[indexOfPoint];
                nFarCells = nFarCellsl[indexOfPoint];
            }

            __threadfence();
            __syncthreads();

            register double gm, gmlin;


            int lock = k / 32 * 32;
            int indexOfFirstPointInWarp = MmortonCodesIdxPointsd[lock];
            int shift = ClosePrefixSuml[indexOfFirstPointInWarp];


            __shared__ double freeSheetIntensity[32];
            __shared__ double freeSheetIntensityLin[32];

            int nBlocks = (nClosePnls + 31) / 32;

            for (int blk = 0; blk < nBlocks; ++blk)
            {
                if (blk * 32 + threadIdx.x < nClosePnls)
                {
                    int chd = closeCellsPfl[shift + (blk * 32 + threadIdx.x)];
                    n = chd - nbodiesd;
                    MmortonCodesIdxd_n = MmortonCodesIdxd[n];

                    freeSheetIntensity[threadIdx.x] = dev_pnl[12 * MmortonCodesIdxd_n + 6];//dev_ptr_freeVortexSheet[MmortonCodesIdxd_n];
                    if (calcLin)
                        freeSheetIntensityLin[threadIdx.x] = dev_pnl[12 * MmortonCodesIdxd_n + 9];//dev_ptr_freeVortexSheetLin[MmortonCodesIdxd_n];
                }

                __syncthreads();

                if (k < npointsd)
                {
                    for (int q = 0; q < 32; ++q)
                    {
                        int pcntr = blk * 32 + q;
                        if (pcntr < nClosePnls)
                        {
                            int chd = closeCellsPfl[shift + pcntr];
                            n = chd - nbodiesd;

                            MmortonCodesIdxd_n = MmortonCodesIdxd[n];

                            double4 begendbuf = *(double4*)(rpnl + MmortonCodesIdxd_n * 4);
                            double2 begbuf = make_double2(begendbuf.x, begendbuf.y);
                            double2 endbuf = make_double2(begendbuf.z, begendbuf.w);

                            gm = freeSheetIntensity[q]; //from sharedMem

                            double tempVelNew;

                            if (indexOfPoint != MmortonCodesIdxd_n)
                            {
                                int position = ClosePrefixSuml[indexOfPoint] + pcntr;

                                double2 i00;
                                i00 = i00save[position];

                                tempVelNew = -gm * (i00.x * tau.x + i00.y * tau.y);
                                val -= tempVelNew;

                                if (vellind)
                                {
                                    gmlin = freeSheetIntensityLin[q]; //from sharedMem
                                    double2 i01, i10, i11;
                                    i01 = i01save[position];
                                    i10 = i10save[position];
                                    i11 = i11save[position];

                                    double tempVelNewB = -gmlin * (i01.x * tau.x + i01.y * tau.y);
                                    val -= tempVelNewB;


                                    double tempVelNewC = -gm * (i10.x * tau.x + i10.y * tau.y);
                                    double tempVelNewD = -gmlin * (i11.x * tau.x + i11.y * tau.y);
                                    vallin -= (tempVelNewC + tempVelNewD);
                                }
                            }//if (posI != posJ)
                        }
                    }
                }//if k

                __syncthreads();
            }
            //*/



            __shared__ double2 momSh[32 * orderAlignment];

            shift = FarPrefixSuml[indexOfFirstPointInWarp];

            nBlocks = (nFarCells + 31) / 32;

            for (int blk = 0; blk < nBlocks; ++blk)
            {
#pragma unroll
                for (int cluster = 0; cluster < 32; ++cluster)
                {
                    if (blk * 32 + cluster < nFarCells)
                    {
                        int chd = farCellsPfl[shift + (blk * 32 + cluster)];
                        register const int srtT = MindexSortTd[chd];
                        mom = momsd + (srtT * orderAlignment);
                        if (diff < order)
                            momSh[cluster * (orderAlignment)+diff] = mom[diff];
                    }
                }

                //*/
                __syncthreads();

                if (k < npointsd)
                {
#pragma unroll
                    for (int qq = 0; qq < 32; ++qq)
                    {
                        int pcntr = blk * 32 + qq;
                        if (pcntr < nFarCells)
                        {
                            int chd = farCellsPfl[shift + pcntr];
                            ps = Mposd[chd];
                            double dist2 = (p.x - ps.x) * (p.x - ps.x) + (p.y - ps.y) * (p.y - ps.y);
                            double2 theta = (p - ps) / dist2;

#pragma unroll
                            for (int q = 0; q < order; ++q)
                                Eloc[q] = make_double2(0.0, 0.0);

#pragma unroll
                            for (int q = 0; q < order; ++q)
                            {
                                for (int s = q; s >= 0; --s)
                                    Eloc[s] += (((2 * q - s) & 1) ? -ifacGPU[q - s + 1] : ifacGPU[q - s + 1]) * multzA(theta, momSh[qq * (orderAlignment)+q - s]);
                                theta = ((q + 1) / dist2) * multz(theta, p - ps);
                            }

                            double2 v = 0.5 * Eloc[0];
                            double2 vL{ 0.0, 0.0 };

                            double2 rPan = di;

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
                                    v += ifacGPU[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                else
                                    if (calcLin)
                                        vL += ifacGPU[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
                            }
                            v *= 2.0;
                            if (calcLin)
                                vL *= 2.0;

                            val += (-v.y * rPan.x + v.x * rPan.y);

                            if (calcLin)
                                vallin += (-vL.y * rPan.x + vL.x * rPan.y);
                        }
                    }//for blk
                }//if k

                __syncthreads();
            }

            // update velocity
            if (k < npointsd)
            {
                double cft = idpid * idlen;
                veld[indexOfPoint] = cft * val;

                if (vellind)
                    vellind[indexOfPoint] = cft * vallin;
            }
        }
    }//treeMatrToVecNonZeroIterKernel12(...)


}





	
#endif
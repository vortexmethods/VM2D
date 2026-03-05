
    inline void Summ15()
    {
        using double4 = Point4D;
        using double2 = Point2D;
        using int2 = std::pair<int, int>;

        //цикл по внутренним узлам снизу вверх

#pragma omp parallel
        {
            int i, j, ch;

            double4 lu[2]; //bounding boxes двух детей

            double2 mom0;
            double2 mom1;
            double2 mom2;
            double2 mom3;
            double2 mom4;
            double2 mom5;
            double2 mom6;
            double2 mom7;
            double2 mom8;
            double2 mom9;
            double2 mom10;
            double2 mom11;
            double2 mom12;
            double2 mom13;
            double2 mom14;
            double2 cen; //центр текущего родительского узла
            double2 dr;  //вектор от центра родителя к центру ребенка

            int m[2]; //для листа это 1, для внутреннего узла это число листьев в нем
            int cm; //масса текущего узла = сумма масс двух детей
            const int nnodes = 2 * (int)object.size() - 1;
            const int nbodies = (int)object.size();
#pragma omp for schedule(dynamic,5)
            for (int k = nbodies; k < nnodes; ++k)
            {
                //MortonTree:
                // 0 1 2 ... (nb-2) x (nb+0) (nb+1) (nb+2) ... (nb+(nb-1))
                // ----------------   -----------------------------------
                //      cells                         bodies

                //Martin's tree:
                // 0 1 2 ... (nb-1) x x x x (nn-(nb-1)) ... (nn-2) (nn-1)
                // ----------------          ----------------------------
                //      bodies                 sorted and reversed cells
                //printf("k = %d\n", k);

                j = 0;
                cm = 0;

                // iterate over all cells assigned to thread
                while (cm == 0)
                {
                    j = 2;
                    int srt = indexSort[(nnodes - 1) - k]; //проход снизу вверх - т.е. в обратном порядке
                    int2 chdPair = child[srt];

                    for (i = 0; i < 2; i++) {
                        int chd = i * chdPair.second + (1 - i) * chdPair.first;   // i==0 => .x;  i==1 => .y

                        ch = (chd >= nbodies) ? (chd - nbodies) : ((nnodes - 1) - indexSortT[chd]);

                        if ((chd >= nbodies) || (mass[nnodes - 1 - ch] >= 0))
                            j--;
                    }

                    if (j == 0)
                    {
                        // all children are ready
                        const int kch = ((nnodes - 1) - k) * orderAlignment; //позиция, куда будут записаны мм для текущего внутреннего узла в массиве moms
                        //moms хранится не по индексу srt, а по отсортированному порядку узлов

    //Считаем bounding box текущего узла из двух детей
                        for (i = 0; i < 2; i++)
                        {
                            const int chd = i * chdPair.second + (1 - i) * chdPair.first; // индекс i-го ребенка
                            if (chd >= nbodies)//если ребенок - это лист
                            {
                                ch = chd - nbodies; //номер частицы в отсортированном массиве
                                const int sortedBody = mortonCodesIdx[ch];

                                double4 xyAB = gabForLeaves[sortedBody];
                                lu[i] = double4{
                                    ::fmin(xyAB[0], xyAB[2]), ::fmin(xyAB[1], xyAB[3]),
                                    ::fmax(xyAB[0], xyAB[2]), ::fmax(xyAB[1], xyAB[3])
                                };
                            }
                            else
                            {
                                ch = indexSortT[chd];
                                // если внутренний узел, его bounding box уже должен быть посчитан, тк идем снизу вверх
                                lu[i] = lowerupper[chd];
                            }
                        }//for i

                        // Объединяем bounding box двух детей:
                        lowerupper[srt] = double4{
                            ::fmin(lu[0][0], lu[1][0]), ::fmin(lu[0][1], lu[1][1]),
                            ::fmax(lu[0][2], lu[1][2]), ::fmax(lu[0][3], lu[1][3])
                        };
                        // Центр текущего узла
                        cen = center[srt] = double2{ 0.5 * (lowerupper[srt][0] + lowerupper[srt][2]),
                                                     0.5 * (lowerupper[srt][1] + lowerupper[srt][3]) };

                        const double2 zero = { 0.0, 0.0 };
                        double2 momh0 = zero;
                        double2 momh1 = zero;
                        double2 momh2 = zero;
                        double2 momh3 = zero;
                        double2 momh4 = zero;
                        double2 momh5 = zero;
                        double2 momh6 = zero;
                        double2 momh7 = zero;
                        double2 momh8 = zero;
                        double2 momh9 = zero;
                        double2 momh10 = zero;
                        double2 momh11 = zero;
                        double2 momh12 = zero;
                        double2 momh13 = zero;
                        double2 momh14 = zero;
                        //переносим мультипольные моменты в центр родительской ячейки
                        for (i = 0; i < 2; i++)
                        {
                            const int chd = i * chdPair.second + (1 - i) * chdPair.first;
                            if (chd >= nbodies) //если ребенок - это лист
                            {
                                ch = chd - nbodies;
                                const int sortedBody = mortonCodesIdx[ch];

                                if (objectType == object_T::point4)
                                {
                                    mom0 = double2{ gamma[sortedBody], 0.0 };
                                    //для вихря все остальные мм нулевые
mom1 = mom2 = mom3 = mom4 = mom5 = mom6 = mom7 = mom8 = mom9 = mom10 = mom11 = mom12 = mom13 = mom14 = double2{ 0.0, 0.0 };
                                    double2 pos = object[sortedBody];
                                    dr = pos - cen;
                                    m[i] = 1;
                                } //objectType==point4

                                if (objectType == object_T::panel)
								{
								  //...
                                } //objectType == panel

                            }
                            else // если ребенок - внутренний узел
                            {
                                const int srtT = indexSortT[chd];
                                ch = (nnodes - 1) - srtT;

                                const int nch = srtT * orderAlignment;

                                mom0 = double2{ moms[nch + 0][0], (double)0 };
                                mom1 = moms[nch + 1];
                                mom2 = moms[nch + 2];
                                mom3 = moms[nch + 3];
                                mom4 = moms[nch + 4];
                                mom5 = moms[nch + 5];
                                mom6 = moms[nch + 6];
                                mom7 = moms[nch + 7];
                                mom8 = moms[nch + 8];
                                mom9 = moms[nch + 9];
                                mom10 = moms[nch + 10];
                                mom11 = moms[nch + 11];
                                mom12 = moms[nch + 12];
                                mom13 = moms[nch + 13];
                                mom14 = moms[nch + 14];
                                dr = center[chd] - cen;
                                m[i] = mass[srtT];
                            }
                            // add child's contribution

                            momh0 += mom0;
                            momh1 += mom1;
                            momh2 += mom2;
                            momh3 += mom3;
                            momh4 += mom4;
                            momh5 += mom5;
                            momh6 += mom6;
                            momh7 += mom7;
                            momh8 += mom8;
                            momh9 += mom9;
                            momh10 += mom10;
                            momh11 += mom11;
                            momh12 += mom12;
                            momh13 += mom13;
                            momh14 += mom14;

                            double2 z = dr;

                            momh1 += multz(mom0, z);
                            momh2 += 2 * multz(mom1, z);
                            momh3 += 3 * multz(mom2, z);
                            momh4 += 4 * multz(mom3, z);
                            momh5 += 5 * multz(mom4, z);
                            momh6 += 6 * multz(mom5, z);
                            momh7 += 7 * multz(mom6, z);
                            momh8 += 8 * multz(mom7, z);
                            momh9 += 9 * multz(mom8, z);
                            momh10 += 10 * multz(mom9, z);
                            momh11 += 11 * multz(mom10, z);
                            momh12 += 12 * multz(mom11, z);
                            momh13 += 13 * multz(mom12, z);
                            momh14 += 14 * multz(mom13, z);

                            z = multz(z, dr);

                            momh2 += multz(mom0, z);
                            momh3 += 3 * multz(mom1, z);
                            momh4 += 6 * multz(mom2, z);
                            momh5 += 10 * multz(mom3, z);
                            momh6 += 15 * multz(mom4, z);
                            momh7 += 21 * multz(mom5, z);
                            momh8 += 28 * multz(mom6, z);
                            momh9 += 36 * multz(mom7, z);
                            momh10 += 45 * multz(mom8, z);
                            momh11 += 55 * multz(mom9, z);
                            momh12 += 66 * multz(mom10, z);
                            momh13 += 78 * multz(mom11, z);
                            momh14 += 91 * multz(mom12, z);

                            z = multz(z, dr);

                            momh3 += multz(mom0, z);
                            momh4 += 4 * multz(mom1, z);
                            momh5 += 10 * multz(mom2, z);
                            momh6 += 20 * multz(mom3, z);
                            momh7 += 35 * multz(mom4, z);
                            momh8 += 56 * multz(mom5, z);
                            momh9 += 84 * multz(mom6, z);
                            momh10 += 120 * multz(mom7, z);
                            momh11 += 165 * multz(mom8, z);
                            momh12 += 220 * multz(mom9, z);
                            momh13 += 286 * multz(mom10, z);
                            momh14 += 364 * multz(mom11, z);

                            z = multz(z, dr);

                            momh4 += multz(mom0, z);
                            momh5 += 5 * multz(mom1, z);
                            momh6 += 15 * multz(mom2, z);
                            momh7 += 35 * multz(mom3, z);
                            momh8 += 70 * multz(mom4, z);
                            momh9 += 126 * multz(mom5, z);
                            momh10 += 210 * multz(mom6, z);
                            momh11 += 330 * multz(mom7, z);
                            momh12 += 495 * multz(mom8, z);
                            momh13 += 715 * multz(mom9, z);
                            momh14 += 1001 * multz(mom10, z);

                            z = multz(z, dr);

                            momh5 += multz(mom0, z);
                            momh6 += 6 * multz(mom1, z);
                            momh7 += 21 * multz(mom2, z);
                            momh8 += 56 * multz(mom3, z);
                            momh9 += 126 * multz(mom4, z);
                            momh10 += 252 * multz(mom5, z);
                            momh11 += 462 * multz(mom6, z);
                            momh12 += 792 * multz(mom7, z);
                            momh13 += 1287 * multz(mom8, z);
                            momh14 += 2002 * multz(mom9, z);

                            z = multz(z, dr);

                            momh6 += multz(mom0, z);
                            momh7 += 7 * multz(mom1, z);
                            momh8 += 28 * multz(mom2, z);
                            momh9 += 84 * multz(mom3, z);
                            momh10 += 210 * multz(mom4, z);
                            momh11 += 462 * multz(mom5, z);
                            momh12 += 924 * multz(mom6, z);
                            momh13 += 1716 * multz(mom7, z);
                            momh14 += 3003 * multz(mom8, z);

                            z = multz(z, dr);

                            momh7 += multz(mom0, z);
                            momh8 += 8 * multz(mom1, z);
                            momh9 += 36 * multz(mom2, z);
                            momh10 += 120 * multz(mom3, z);
                            momh11 += 330 * multz(mom4, z);
                            momh12 += 792 * multz(mom5, z);
                            momh13 += 1716 * multz(mom6, z);
                            momh14 += 3432 * multz(mom7, z);

                            z = multz(z, dr);

                            momh8 += multz(mom0, z);
                            momh9 += 9 * multz(mom1, z);
                            momh10 += 45 * multz(mom2, z);
                            momh11 += 165 * multz(mom3, z);
                            momh12 += 495 * multz(mom4, z);
                            momh13 += 1287 * multz(mom5, z);
                            momh14 += 3003 * multz(mom6, z);

                            z = multz(z, dr);

                            momh9 += multz(mom0, z);
                            momh10 += 10 * multz(mom1, z);
                            momh11 += 55 * multz(mom2, z);
                            momh12 += 220 * multz(mom3, z);
                            momh13 += 715 * multz(mom4, z);
                            momh14 += 2002 * multz(mom5, z);

                            z = multz(z, dr);

                            momh10 += multz(mom0, z);
                            momh11 += 11 * multz(mom1, z);
                            momh12 += 66 * multz(mom2, z);
                            momh13 += 286 * multz(mom3, z);
                            momh14 += 1001 * multz(mom4, z);

                            z = multz(z, dr);

                            momh11 += multz(mom0, z);
                            momh12 += 12 * multz(mom1, z);
                            momh13 += 78 * multz(mom2, z);
                            momh14 += 364 * multz(mom3, z);

                            z = multz(z, dr);

                            momh12 += multz(mom0, z);
                            momh13 += 13 * multz(mom1, z);
                            momh14 += 91 * multz(mom2, z);

                            z = multz(z, dr);

                            momh13 += multz(mom0, z);
                            momh14 += 14 * multz(mom1, z);

                            z = multz(z, dr);

                            momh14 += multz(mom0, z);
                        }

                        // Сохраняем итоговые моменты текущего внутреннего узла в общий массив
                        momh0[1] = 0;
                        moms[kch + 0] = momh0;
                        moms[kch + 1] = momh1;
                        moms[kch + 2] = momh2;
                        moms[kch + 3] = momh3;
                        moms[kch + 4] = momh4;
                        moms[kch + 5] = momh5;
                        moms[kch + 6] = momh6;
                        moms[kch + 7] = momh7;
                        moms[kch + 8] = momh8;
                        moms[kch + 9] = momh9;
                        moms[kch + 10] = momh10;
                        moms[kch + 11] = momh11;
                        moms[kch + 12] = momh12;
                        moms[kch + 13] = momh13;
                        moms[kch + 14] = momh14;
                        cm = m[0] + m[1];
                    }

#pragma omp flush

                    if (cm != 0)
                    {
                        mass[nnodes - 1 - k] = cm;// Записываем массу текущего узла
                        //k += inc;
                        //flag = 0;
                    }
                }//while flag==0
            }//for k
        }
    }//Summ12()
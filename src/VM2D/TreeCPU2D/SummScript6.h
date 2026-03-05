
    inline void Summ6()
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
mom1 = mom2 = mom3 = mom4 = mom5 = double2{ 0.0, 0.0 };
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

                            double2 z = dr;

                            momh1 += multz(mom0, z);
                            momh2 += 2 * multz(mom1, z);
                            momh3 += 3 * multz(mom2, z);
                            momh4 += 4 * multz(mom3, z);
                            momh5 += 5 * multz(mom4, z);

                            z = multz(z, dr);

                            momh2 += multz(mom0, z);
                            momh3 += 3 * multz(mom1, z);
                            momh4 += 6 * multz(mom2, z);
                            momh5 += 10 * multz(mom3, z);

                            z = multz(z, dr);

                            momh3 += multz(mom0, z);
                            momh4 += 4 * multz(mom1, z);
                            momh5 += 10 * multz(mom2, z);

                            z = multz(z, dr);

                            momh4 += multz(mom0, z);
                            momh5 += 5 * multz(mom1, z);

                            z = multz(z, dr);

                            momh5 += multz(mom0, z);
                        }

                        // Сохраняем итоговые моменты текущего внутреннего узла в общий массив
                        momh0[1] = 0;
                        moms[kch + 0] = momh0;
                        moms[kch + 1] = momh1;
                        moms[kch + 2] = momh2;
                        moms[kch + 3] = momh3;
                        moms[kch + 4] = momh4;
                        moms[kch + 5] = momh5;
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
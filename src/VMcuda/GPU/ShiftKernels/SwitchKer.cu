switch (order)
{
case 1:
    SummarizationKernel2_1<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 2:
    SummarizationKernel2_2<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 3:
    SummarizationKernel2_3<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 4:
    SummarizationKernel2_4<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 5:
    SummarizationKernel2_5<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 6:
    SummarizationKernel2_6<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 7:
    SummarizationKernel2_7<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 8:
    SummarizationKernel2_8<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 9:
    SummarizationKernel2_9<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 10:
    SummarizationKernel2_10<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 11:
    SummarizationKernel2_11<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 12:
    SummarizationKernel2_12<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 13:
    SummarizationKernel2_13<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 14:
    SummarizationKernel2_14<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 15:
    SummarizationKernel2_15<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

case 16:
    SummarizationKernel2_16<<<blocks * FACTOR3, THREADS3>>>(nnodesd, nbodiesd, (int2*)Mchildd, massd, order, (real2*)momsd, (double*)vtxd, objType, MmortonCodesIdxd, (real2*)Mposd, MindexSortd, MindexSortTd);
    break;

}
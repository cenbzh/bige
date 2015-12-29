/*************************************************************************
	> File Name: bige_output_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 19时53分30秒
 ************************************************************************/

#include<stdio.h>
#include "bige_output_module.h"
#include "bige_external.h"
void bige_output(Population* pop,FILE* fp)
{
    int i,j;
    Individual* ind;
    int count=popsize;
    int printflag=1;
    for(i=0;i<popsize;i++)
    {
        if(printflag==1)
        {
            ind=&(pop->ind[i]);
            for(j=0;j<nfunc;j++)
            {
                fprintf(fp,"%12lf",ind->objs[j]);
            }
            fprintf(fp,"\n");
        }
    }
    return;
}

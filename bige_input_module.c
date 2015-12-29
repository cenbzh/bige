/*************************************************************************
	> File Name: bige_input_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 19时40分55秒
 ************************************************************************/

#include<stdio.h>
#include "bige_input_module.h"
#include "bige_external.h"


void bige_input(FILE* fp)
{
    /*int generNum;
    int varNum;
    int funcNum;
    double crossp;
    double mutp;
    double lim_b;
    double lim_u;
    int popnum;
    double ditemp;
    double dimtemp;
    int ans;*/
    int i;
    fscanf(fp,"%d",&ngener);/*get no of generations*/
    fscanf(fp,"%d",&popsize);
    fscanf(fp,"%d",&nvar);
    fscanf(fp,"%d",&nfunc);
    fscanf(fp,"%lf",&pcross);
    fscanf(fp,"%lf",&pmut_r);
    fscanf(fp,"%d",&isLimit);
    if(isLimit==1)
    {
        for(i=0;i<nvar;i++)
        {
            lim_r[i][0]=0;
            lim_r[i][1]=1;
        }
    }
    fscanf(fp,"%lf",&di);
    fscanf(fp,"%lf",&dim);
    return;
}

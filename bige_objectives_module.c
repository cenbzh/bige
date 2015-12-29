/*************************************************************************
	> File Name: bige_objectives_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 19时12分24秒
 ************************************************************************/

#include<stdio.h>
#include "bige_objectives_module.h"
#include <string.h>
#include "bige_external.h"

#define TESTNO1 1
#define TESTNAME1 "WFG"


extern void wfg_eval(double*,int,int,int,char*,double*);

void bige_getobjects(Population* oldpop_ptr, int size, char* problem,char* testdata)
{
    int testno;
    bige_testdata(testdata,&testno);
    int p,j;
    Individual* ind;
    for(p=0;p<size;p++)
    {
        ind=&(oldpop_ptr->ind[p]);
        switch(testno)
        {
            case TESTNO1:
                wfg_eval(ind->xreal,nvar,2*(nfunc-1),nfunc,problem,ind->objs);
                break;
            default:
                break;
        }
        for(j=0;j<nfunc;j++)
        {
            if(p==0)
            {
                oldpop_ptr->maxObj[j]=ind->objs[j];
                oldpop_ptr->minObj[j]=ind->objs[j];
            }
            else
            {
                if(oldpop_ptr->maxObj[j]<ind->objs[j])
                {
                    oldpop_ptr->maxObj[j]=ind->objs[j];
                }
                if(oldpop_ptr->minObj[j]>=ind->objs[j])
                {
                    oldpop_ptr->minObj[j]=ind->objs[j];
                }
            }
        }
    }
}

void bige_testdata(char* testdata,int* testno)
{
    if(strcmp(testdata,TESTNAME1)==0)
    {
        *testno=TESTNO1;
    }
    return;
}

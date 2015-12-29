/*************************************************************************
	> File Name: bige_main.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 20时12分53秒
 ************************************************************************/

#include<stdio.h>
#include "bige_define.h"
#include "bige_struct.h"
#include "bige_estimation_module.h"
#include "bige_input_module.h"
#include "bige_output_module.h"
#include "bige_objectives_module.h"
#include "bige_operator_module.h"
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define DIR_PATH "/home/cenbzh/Documents/paper/bige/c/obj5/"

int nvar;
int ngener;
int isLimit;
double pcross;
double pmut_r;
int nfunc;
int popsize;
double radiu;
double di,dim;
int shMatrix[2*MAXPOP][2*MAXPOP];
int gpopsize;
Population oldPops;
Population matePops;
Population newPops;
LayerList layerlist;
double maxObjs[MAXFUN];

void bige_engin(FILE* fp,char* problem,char* testdata)
{
    int i;
    Population* oldpop_ptr=&oldPops;
    Population* newpop_ptr=&newPops;
    Population* matepop_ptr=&matePops;
    LayerList* list_ptr=&layerlist;
    radiu=pow(1.0/popsize,1.0/nfunc);
    printf("init pop...\n");
    bige_init_pop(oldpop_ptr,popsize);

    printf("bige start...\n");
    for(i=0;i<ngener;i++)
    {
        printf("gener: %d\n",i+1);
        bige_getobjects(oldpop_ptr,popsize,problem,testdata);
        bige_estimation_pr(oldpop_ptr,popsize);
        bige_estimation_cd(oldpop_ptr,popsize);
        bige_mate_select(oldpop_ptr,matepop_ptr);
        bige_generate_offspring(matepop_ptr,newpop_ptr);
        bige_getobjects(newpop_ptr,popsize,problem,testdata);
        bige_init_list(list_ptr);
        bige_env_select(oldpop_ptr,newpop_ptr,matepop_ptr,list_ptr,problem,testdata);
        bige_copy_pop(oldpop_ptr,matepop_ptr);
        bige_clear_list(list_ptr);
    }
    bige_getobjects(oldpop_ptr,popsize,problem,testdata);
    bige_output(oldpop_ptr,fp);
    for(i=0;i<nfunc;i++)
    {
        if(maxObjs[i]<oldpop_ptr->maxObj[i])
        {
            maxObjs[i]=oldpop_ptr->maxObj[i];
        }
    }
    return;
}

int main(int argc,char* argv[])
{
    if(argc!=3)
    {
        printf("usage: %s <problem> <problemtype>",argv[0]);
        exit(-1);
    }
    int i;
    srand(time(NULL));
    FILE* infile;
    FILE* outfile;
    char inName[200];
    strcpy(inName,DIR_PATH);
    strcat(inName,"input.txt");
    char outName[200];
    char format[10];
    memset(maxObjs,0,sizeof(maxObjs));
    int runtime=1;
    strcpy(outName,DIR_PATH);
    strcat(outName,argv[1]);
    strcat(outName,"finalfit.txt");
    outfile=fopen(outName,"w");
    for(i=0;i<runtime;i++)
    {
        infile=fopen(inName,"r");
        bige_input(infile);
        bige_engin(outfile,argv[1],argv[2]);
        if(i==runtime-1)
        {
            int j;
            for(j=0;j<nfunc;j++)
            {
                fprintf(outfile,"%12lf",maxObjs[j]);
            }
        }
        fclose(infile);
    }
    fclose(outfile);
    return 0;
}

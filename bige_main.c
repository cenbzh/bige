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
#define DIR_PATH "/home/cenbzh/Documents/paper/bige/c/obj10/"

int nvar;
int ngener;
int isLimit;
double lim_r[MAXVAR][2];
double pcross;
double pmut_r;
int nfunc;
int popsize;
double radiu;
double di,dim;
double shMatrix[2*MAXPOP][2*MAXPOP];
int gpopsize;
Population oldPops;
Population matePops;
Population newPops;
LayerList layerlist;
double maxObjs[MAXFUN];


/*用来测试的变量*/
int calln;
int equaln;

void bige_engin(FILE* fp,FILE* fp1,char* problem,char* testdata)
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
        bige_prcd_out(oldpop_ptr,fp1);
        bige_mate_select(oldpop_ptr,matepop_ptr);
        bige_generate_offspring(matepop_ptr,newpop_ptr);
        bige_getobjects(newpop_ptr,popsize,problem,testdata);
        bige_init_list(list_ptr);
        bige_env_select(oldpop_ptr,newpop_ptr,matepop_ptr,list_ptr,problem,testdata);
        bige_copy_pop(oldpop_ptr,matepop_ptr);
        bige_clear_list(list_ptr);
    }
    bige_getobjects(oldpop_ptr,popsize,problem,testdata);
    bige_estimation_pr(oldpop_ptr,popsize);
    bige_estimation_cd(oldpop_ptr,popsize);
    bige_output(oldpop_ptr,fp);
    bige_prcd_out(oldpop_ptr,fp1);
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
    equaln=0;
    calln=0;

    int i;
    srand(time(NULL));
    FILE* infile;
    FILE* outfile;
    FILE* prcdfile;

    char inName[200];
    strcpy(inName,DIR_PATH);
    strcat(inName,"input.txt");
    char outName[200];
    char prcdName[200];
    memset(maxObjs,0,sizeof(maxObjs));
    int runtime=30;
    strcpy(outName,DIR_PATH);
    strcat(outName,argv[1]);
    strcat(outName,"/chebyshev_finalfit1.txt");
    outfile=fopen(outName,"w");

    strcpy(prcdName,DIR_PATH);
    strcat(prcdName,argv[1]);
    strcat(prcdName,"/chebyshev_prcd1.txt");
    prcdfile=fopen(prcdName,"w");
    for(i=0;i<runtime;i++)
    {
        infile=fopen(inName,"r");
        bige_input(infile);
        bige_engin(outfile,prcdfile,argv[1],argv[2]);
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
    fclose(prcdfile);
    fclose(outfile);
    printf("call times: %d, equal times: %d\n",calln, equaln);
    return 0;
}

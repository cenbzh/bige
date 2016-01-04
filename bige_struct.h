/*************************************************************************
	> File Name: bige_struct.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月15日 星期二 23时14分19秒
 ************************************************************************/

#ifndef _BIGE_STRUCT_H
#define _BIGE_STRUCT_H
#include "bige_define.h"

typedef struct Individual
{
    int rank;
    int rank2;
    int dominated;
    double xreal[MAXVAR];
    double objs[MAXFUN];
    double proximity;
    double crowdingDegree;
}Individual;

typedef struct Population
{
    int maxrank;
    double maxObj[MAXFUN];
    double minObj[MAXFUN];
    Individual ind[2*MAXPOP];
}Population;

typedef struct LayerIndSet
{
    int indNo;
    struct LayerIndSet* next;
}LayerIndSet;

typedef struct DominIndSet
{
    int indNo;
    struct DominIndSet* next;
}DominIndSet;

typedef struct LayerList
{
    int nlayer;
    LayerIndSet* layers[2*MAXPOP];
    DominIndSet* dominset[2*MAXPOP];
    int layerIndNum[2*MAXPOP];
}LayerList;

#endif

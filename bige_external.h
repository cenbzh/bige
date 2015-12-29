/*************************************************************************
	> File Name: bige_external.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 17时00分36秒
 ************************************************************************/

#ifndef _BIGE_EXTERNAL_H
#define _BIGE_EXTERNAL_H
#include "bige_define.h"
extern int nfunc;
extern int nvar;
extern int popsize;
extern double pcross;
extern int isLimit;
extern double pmut_r;
extern double radiu;
extern double lim_r[MAXVAR][2];
extern double di,dim;
extern int ngener;
extern double shMatrix[2*MAXPOP][2*MAXPOP];
extern int gpopsize;
#endif

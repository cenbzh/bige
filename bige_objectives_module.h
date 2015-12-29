/*************************************************************************
	> File Name: bige_objectives_module.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 19时08分29秒
 ************************************************************************/

#ifndef _BIGE_OBJECTIVES_MODULE_H
#define _BIGE_OBJECTIVES_MODULE_H
#include "bige_struct.h"
void bige_getobjects(Population* oldpop_ptr,int size,char* problem,char* testdata);
void bige_testdata(char* testdata,int* tno);
#endif

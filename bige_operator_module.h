/*************************************************************************
	> File Name: bige_operator_module.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月16日 星期三 23时42分56秒
 ************************************************************************/

#ifndef _BIGE_OPERATOR_MODULE_H
#define _BIGE_OPERATOR_MODULE_H
#include "bige_struct.h"

void bige_init_pop(Population* pop,int size);/*初始化种群*/

void bige_generate_offspring(Population* oldpop_ptr,Population* newpop_ptr);/*产生后代*/

void bige_mate_select(Population* oldpop_ptr,Population* matepop_ptr);/*繁殖选择*/

void bige_env_select(Population* oldpop_ptr,Population* newpop_ptr,Population* nextpop_ptr, LayerList* list, char* problem);/*环境选择*/

Individual* bige_tournament(Individual* d1,Individual* d2);/*锦标赛选择*/

void bige_copy_pop(Population* des_ptr,Population* scr_ptr);/*复制种群*/

int bige_comp_ind(Individual* d1,Individual* d2);/*判断两个个体的占优关系*/

void bige_copy_ind(Individual* des_ptr,Individual* scr_ptr);/*复制个体*/

void bige_sbx(Population* oldpop_ptr,Population* newpop_ptr);/*交叉操作*/

void bige_mutation(Population* newpop_ptr);/*突变操作*/

void bige_union_pop(Population* pop1,Population* pop2);/*合并种群*/

void bige_nondominate_sort(Population* pop,LayerList* list);

void bige_keepalive(Population* pop, LayerList* list, Population* nextpop_ptr);

void bige_clear_list(LayerList* list);

void bige_init_list(LayerList* list);

void bige_add_layerind(LayerList* list,int i,int n);

void bige_add_dominind(LayerList* list,int i,int j);

#endif

/*************************************************************************
	> File Name: bige_estimation_module.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月15日 星期二 20时00分48秒
 ************************************************************************/

#ifndef _BIGE_ESTIMATION_MODULE_H
#define _BIGE_ESTIMATION_MODULE_H
#include "bige_struct.h"

void bige_estimation_pr(Population* pop,int size);/*对种群pop的个体进行收敛度的评估*/

void bige_estimation_cd(Population* pop,int size);/*对种群pop的个体进行收敛度的评估*/

double bige_share_ind_cd(Population* pop,int p,int size);/*获得下标为p的个体的拥挤度*/

double bige_distance(Population* pop,int p,int q);/*计算下标为p的个体与下标为q的个体的距离*/

double bige_value_normalization(Population* pop,int p,int k);/*对下标p的个体的第k个目标进行归一化*/

void bige_share_function(Population* pop,int size);/*计算种群个体的共享值*/

double bige_chebyshev(Population* pop,int p);

double bige_sum(Population* pop,int p);

double bige_angle_ind_cd(Population* pop,int p,int size);

double bige_compute_angle(Population* pop,int p,int q);

double bige_angle_assign(Population* pop,int size);

double bige_angle_ind_cd_kclosest(Population* pop,int p,int size);

#endif

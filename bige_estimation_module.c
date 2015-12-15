/*************************************************************************
	> File Name: bige_estimation_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月15日 星期二 20时21分33秒
 ************************************************************************/

#include "bige_estimation_module.h"
#include <math.h>
#include "bige_difine.h"
/*这里可以选择不同的计算距离方法，示例为欧氏距离*/

extern int nfunc;
extern int shMatrix[2*MAXPOP][2*MAXPOP];
double bige_distance(Population* pop,int p,int q)
{

    int k;
    double dis=0;
    for(k=0;k<nfunc;k++)
    {
        double fp=bige_value_normalization(pop,p,k);
        double fq=bige_value_normalization(pop,q,k);
        dis+=(fp-fq)*(fp-fq);
    }
    return pow(dis,1.0/2);
}

/*这里可以选择不同的评估方法，示例为论文的评估方法*/
void bige_estimation_pr(Population* pop,int size)
{
    int p,k;
    double sum=0;
    Individual* ind;
    for(p=0;p<size;p++)//计算种群个体的收敛度
    {
        sum+=bige_value_normalization(pop,p,k);
    }
    ind->proximity=sum;
    return;
}


/*这里可以选择不同的评估方法，示例为论文的评估方法*/
void bige_estimation_cd(Population* pop,int size)
{
    
}
double bige_get_ind_cd(Population* pop,int p,int size)
{
}

/*这里可以选择不同的归一化方法，示例为比较常见的方法*/
double bige_value_normalization(Population* pop,int p,int k)
{
    if(pop->maxObj[k]==pop->minObj[k])
    {
        return 0;
    }
    Individual* ind=&(pop->ind[p]);
    return (ind->objs[k]-pop->minObj[k])/(pop->maxObj[k]-pop->minObj[k]);
}

void bige_share_function(Population* pop,int size)
{

}

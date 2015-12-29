/*************************************************************************
	> File Name: bige_estimation_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月15日 星期二 20时21分33秒
 ************************************************************************/

#include "bige_estimation_module.h"
#include <math.h>
#include "bige_define.h"
#include <string.h>
#include <stdlib.h>
/*这里可以选择不同的计算距离方法，示例为欧氏距离*/

extern int nfunc;
extern double shMatrix[2*MAXPOP][2*MAXPOP]i;
extern double radiu;
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
        sum=0;
        ind=&(pop->ind[p]);
        for(k=0;k<nfunc;k++)
        {
            sum+=bige_value_normalization(pop,p,k);
        }
        ind->proximity=sum;
    }
    return;
}


/*这里可以选择不同的评估方法，示例为论文的评估方法*/
void bige_estimation_cd(Population* pop,int size)
{
    bige_share_function(pop,size);
    int p,k,q;
    Individual* ind;
    for(p=0;p<size;p++)
    {
        ind=&(pop->ind[p]);
        ind->crowdingDegree=bige_get_ind_cd(pop,p,size);
    }
    return;

}
double bige_get_ind_cd(Population* pop,int p,int size)
{
    int q;
    double sum=0;
    for(q-0;q<size;q++)
    {
        if(q==p)
        {
            continue;
        }
        double sh=shMatrix[p][q];
        sum+=sh;
    }
    return pow(sum,1.0/2);
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
    memset(shMatrix,0,sizeof(shMatrix));
    int i,j;
    Individual* indi;
    Individual* indj;
    double dis=0;
    double pri;
    double prj;
    double rnk;
    for(i=0;i<size;i++)
    {
        indi=&(pop->ind[i]);
        for(j=i+1;j<size;j++)
        {
            indj=&(pop->ind[j]);
            dis=bige_distance(pop,i,j);
            if(dis>radiu)
            {
                shMatrix[i][j]=0;
                shMatrix[j][i]=0;
            }
            else
            {
                pri=indi->proximity;
                prj=indj->proximity;
                if(pri==prj)
                {
                    rnk=rand()/(1.0*RAND_MAX);
                    if(rnk==0.5)
                    {
                        shMatrix[i][j]=pow(0.5*(1-dis/radiu),2);
                        shMatrix[j][i]=pow(1.5*(1-dis/radiu),2);
                    }
                    else
                    {
                        shMatrix[i][j]=pow(1.5*(1-dis/radiu),2);
                        shMatrix[j][i]=pow(0.5*(1-dis/radiu),2);
                    }
                }
                else if(pri<prj)
                {
                    shMatrix[i][j]=pow(0.5*(1-dis/radiu),2);
                    shMatrix[j][i]=pow(1.5*(1-dis/radiu),2);
                }
                else
                {
                    shMatrix[j][i]=pow(0.5*(1-dis/radiu),2);
                    shMatrix[i][j]=pow(1.5*(1-dis/radiu),2);
                }
            }
        }
    }
    return;
}

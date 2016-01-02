/*************************************************************************
	> File Name: bige_estimation_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月15日 星期二 20时21分33秒
 ************************************************************************/
#include <stdio.h>
#include "bige_estimation_module.h"
#include <math.h>
#include "bige_define.h"
#include <string.h>
#include <stdlib.h>
#include "bige_external.h"

#define PI 3.1415926
extern int calln;
extern int equaln;
//static double weight[2*MAXPOP][MAXFUN];
static double angles[2*MAXPOP][2*MAXPOP];
static int flags[2*MAXPOP];
/*这里可以选择不同的计算距离方法，示例为欧氏距离*/

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
    Individual* ind;
    for(p=0;p<size;p++)//计算种群个体的收敛度
    {
        ind=&(pop->ind[p]);
        ind->proximity=bige_chebyshev(pop,p);
    }
    return;
}

/*这里可以选择不同的评估方法，示例为论文的评估方法*/
void bige_estimation_cd(Population* pop,int size)
{
    //bige_share_function(pop,size);
    bige_angle_assign(pop,size);
    int p,k,q;
    Individual* ind;
    for(p=0;p<size;p++)
    {
        ind=&(pop->ind[p]);
        //ind->crowdingDegree=bige_share_ind_cd(pop,p,size);
        //ind->crowdingDegree=bige_angle_ind_cd(pop,p,size);
        ind->crowdingDegree=5*PI-bige_angle_ind_cd_kclosest(pop,p,size);
    }
    return;

}
double bige_share_ind_cd(Population* pop,int p,int size)
{
    int q;
    double sum=0;
    for(q=0;q<size;q++)
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


/*使用切比雪夫距离*/
double bige_chebyshev(Population* pop,int p)
{
    //calln++;
    int j;
    //memset(weight,0,sizeof(weight));
    Individual* ind;
    ind=&(pop->ind[p]);
    double sum=0;
    for(j=0;j<nfunc;j++)
    {
        if(ind->objs[j]==pop->minObj[j])
        {
            //equaln++;
            return 0;
        }
        else
        {
            sum+=1.0/(ind->objs[j]-pop->minObj[j]);
        }
    }
    return 1.0/sum;
}

/*求和*/
double bige_sum(Population* pop,int p)
{
    int j;
    double sum=0;
    for(j=0;j<nfunc;j++)
    {
        sum+=bige_value_normalization(pop,p,j);
    }
    return sum;
}

/*计算种群个体的角度*/
double bige_compute_angle(Population* pop,int p,int q)
{
    double sum1=0;
    double sum2=0;
    double sum3=0;
    int j;
    Individual* ind1=&(pop->ind[p]);
    Individual* ind2=&(pop->ind[q]);
    for(j=0;j<nfunc;j++)
    {
        sum1+=(ind1->objs[j]-pop->minObj[j])*(ind2->objs[j]-pop->minObj[j]);
        sum2+=(ind1->objs[j]-pop->minObj[j])*(ind1->objs[j]-pop->minObj[j]);
        sum3+=(ind2->objs[j]-pop->minObj[j])*(ind2->objs[j]-pop->minObj[j]);
    }
    return acos(sum1/pow(sum2*sum3,1.0/2));
}

double bige_angle_ind_cd(Population* pop,int p,int size)
{
    double limit=PI/5;
    int i,j;
    Individual *indp;
    indp=&(pop->ind[p]);
    double count=0;
    double angle;
    for(i=0;i<size;i++)
    {
        if(i==p)
        {
            continue;
        }
        calln++;
        angle=bige_compute_angle(pop,p,i);
        //printf("%lf\n",angle);
        if(angle<=limit)
        {
            count++;
            equaln++;
        }
    }
    return count;
}

/*选出k个最小值加起来,然后作为多样性的评估值*/

double bige_angle_ind_cd_kclosest(Population* pop,int p,int size)
{
    int k;
    double sum=0;
    int i;
    int mini;
    memset(flags,0,sizeof(flags));
    for(k=0;k<10;k++)
    {
        double min=-1;
        for(i=0;i<size;i++)
        {
            if(i==p||flags[i]==1)
            {
                continue;
            }
            if(min==-1)
            {
                min=angles[p][i];
                mini=i;
            }
            else
            {
                if(min>angles[p][i])
                {
                    min=angles[p][i];
                    mini=i;
                }
            }
        }
        sum+=min;
        flags[mini]=1;
    }
    return sum;
}

double bige_angle_assign(Population* pop,int size)
{
    memset(angles,0,sizeof(angles));
    int i,j;
    double angle;
    for(i=0;i<size;i++)
    {
        for(j=i+1;j<size;j++)
        {
            angle=bige_compute_angle(pop,i,j);
            angles[i][j]=angle;
            angles[j][i]=angle;
        }
    }
    return 0;
}

/*************************************************************************
	> File Name: bige_operator_module.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年12月29日 星期二 15时06分50秒
 ************************************************************************/

#include<stdio.h>
#include "bige_operator_module.h"
#include "bige_define.h"
#include "bige_struct.h"
#include "bige_external.h"
#include "bige_estimation_module.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bige_objectives_module.h"

void bige_init_pop(Population* pop,int size)
{
    int i,j;
    double d;
    pop->maxrank=0;
    Individual* ind=&(pop->ind[0]);
    for(i=0;i<popsize;i++)
    {
        ind->rank=0;
        ind->dominated=0;
        ind->proximity=0;
        ind->crowdingDegree=0;
        for(j=0;j<nvar;j++)
        {
            d=rand()/(1.0*RAND_MAX);;
            ind->xreal[j]=d;
        }
        ind++;
    }
    return;
}

void bige_generate_offspring(Population* matepop_ptr,Population* newpop_ptr)
{
    bige_sbx(matepop_ptr,newpop_ptr);
    bige_mutation(newpop_ptr);
    return;
}

void bige_mate_select(Population* oldpop_ptr,Population* matepop_ptr)
{
    int ind1,ind2;
    int p,k;
    Individual* select_ptr;
    Individual* ind1_ptr;
    Individual* ind2_ptr;
    Individual* mateind_ptr;
    Individual* lastselect=NULL;
    for(p=0,k=0;p<popsize;)
    {
        ind1=rand()%popsize;
        ind2=rand()%popsize;
        while(ind2==ind1)
        {
            ind2=rand()%popsize;
        }
        ind1_ptr=&(oldpop_ptr->ind[ind1]);
        ind2_ptr=&(oldpop_ptr->ind[ind2]);
        mateind_ptr=&(matepop_ptr->ind[k]);
        select_ptr=bige_tournament(ind1_ptr,ind2_ptr);
        if(select_ptr==lastselect&& p%2==1)
        {
            continue;
        }
        bige_copy_ind(mateind_ptr,select_ptr);
        lastselect=select_ptr;
        p++;
        k++;
    }
    return;
}

void bige_env_select(Population* oldpop_ptr,Population* newpop_ptr,Population* nextpop_ptr, LayerList* list, char* problem,char* testdata)
{
    bige_union_pop(oldpop_ptr,newpop_ptr);
    bige_getobjects(oldpop_ptr,2*popsize,problem,testdata);
    gpopsize=2*popsize;
    bige_estimation_pr(oldpop_ptr,gpopsize);
    bige_estimation_cd(oldpop_ptr,gpopsize);
    bige_nondominate_sort(oldpop_ptr,list);
    bige_keepalive(oldpop_ptr,list,nextpop_ptr);
    return;
}

Individual* bige_tournament(Individual* d1,Individual* d2)
{
    int r=bige_comp_ind(d1,d2);
    if(r>0)
    {
        return d1;
    }
    else if(r<0)
    {
        return d2;
    }
    else
    {
        double pr=rand()/(1.0*RAND_MAX);
        if(pr<=0.5)
        {
            return d1;
        }
        else
        {
            return d2;
        }
    }
}

void bige_copy_pop(Population* des_ptr,Population* scr_ptr)
{
    int p;
    Individual* ind1;
    Individual* ind2;
    des_ptr->maxrank=scr_ptr->maxrank;
    for(p=0;p<popsize;p++)
    {
        ind1=&(des_ptr->ind[p]);
        ind2=&(scr_ptr->ind[p]);
        bige_copy_ind(ind1,ind2);
    }
    int i;
    for(i=0;i<nfunc;i++)
    {
        des_ptr->maxObj[i]=scr_ptr->maxObj[i];
        des_ptr->minObj[i]=scr_ptr->minObj[i];
    }
    return ;
}

int bige_comp_ind(Individual* ind1_ptr,Individual* ind2_ptr)
{
    double pr1,pr2,cd1,cd2;
    pr1=ind1_ptr->proximity;
    pr2=ind2_ptr->proximity;
    cd1=ind1_ptr->crowdingDegree;
    cd2=ind2_ptr->crowdingDegree;
    if((pr1<pr2 && cd1<=cd2) || (pr1 <= pr2 && cd1 < cd2))
    {
        return 1;
    }
    else if((pr1>pr2 && cd1>=cd2) || (pr1>=pr2 && cd1 > cd2))
    {
        return -1;
    }
    return 0;
}

void bige_copy_ind(Individual* des_ptr,Individual* scr_ptr)
{
    des_ptr->rank=scr_ptr->rank;
    des_ptr->dominated=scr_ptr->dominated;
    int i;
    for(i=0;i<nvar;i++)
    {
        des_ptr->xreal[i]=scr_ptr->xreal[i];
    }
    for(i=0;i<nfunc;i++)
    {
        des_ptr->objs[i]=scr_ptr->objs[i];
    }
    des_ptr->proximity=scr_ptr->proximity;
    des_ptr->crowdingDegree=scr_ptr->crowdingDegree;
    return;
}

void bige_sbx(Population* matepop_ptr, Population* newpop_ptr)
{
    int i,j,y,n,r;
    double rnd,par1,par2,chld1,chld2,betaq,beta,alpha;
    double y1,y2,yu,yl,expp;
    y=0; n=0;
    for(i = 0; i < popsize/2; i++)
    {
        rnd = rand()/(1.0*RAND_MAX);
        //rnd=randomperc();
        if(rnd <= pcross)
        {
            for(j = 0;j < nvar;j++)
            { 
                par1 = matepop_ptr->ind[y].xreal[j];
                par2 = matepop_ptr->ind[y+1].xreal[j]; 
                yl = lim_r[j][0];
                yu = lim_r[j][1];
                rnd = rand()/(1.0*RAND_MAX);
                //rnd=randomperc();
                if(rnd <= 0.5)
                {
                    if(fabs(par1 - par2) > 0.000001) // changed by Deb (31/10/01)
                    {
                        if(par2 > par1)
                        {
                            y2 = par2;
                            y1 = par1;
                        }
                        else
                        {
                            y2 = par1;
                            y1 = par2;
                        }
                        if((y1 - yl) > (yu - y2))
                        {
                            beta = 1 + (2*(yu - y2)/(y2 - y1));
         	                //printf("beta = %f\n",beta);
                        }
                        else
                        {
                            beta = 1 + (2*(y1-yl)/(y2-y1));
                            //printf("beta = %f\n",beta);
                        }
                        expp = di + 1.0;
                        beta = 1.0/beta;
                        alpha = 2.0 - pow(beta,expp);
                        if (alpha < 0.0) 
                        {
                            printf("ERRRROR %f %d %d %f %f\n",alpha,y,n,par1,par2);
                            exit(-1);
                        }
                        rnd = rand()/(1.0*RAND_MAX); 
                        //rnd=randomperc();
                        if (rnd <= 1.0/alpha)
                        {
                            alpha = alpha*rnd;
                            expp = 1.0/(di+1.0);
                            betaq = pow(alpha,expp);
                        }
                        else
                        {
                            alpha = alpha*rnd;
                            alpha = 1.0/(2.0-alpha);
                            expp = 1.0/(di+1.0);
                            if (alpha < 0.0) 
                            {
                                printf("ERRRORRR \n");
                                exit(-1);
                            }
                            betaq = pow(alpha,expp);
                        }
                        chld1 = 0.5*((y1+y2) - betaq*(y2-y1));
                        chld2 = 0.5*((y1+y2) + betaq*(y2-y1));
                    }
                    else
                    {
                        betaq = 1.0;
                        y1 = par1; y2 = par2;
                        chld1 = 0.5*((y1+y2) - betaq*(y2-y1));
                        chld2 =  0.5*((y1+y2) + betaq*(y2-y1));
                    }
                    //if (chld1 < yl) chld1 = yl;
                    /*while(chld1<yl)
                    {
                        chld1+=(yu-yl);
                    }
                    //if (chld1 > yu) chld1 = yu;
                    while(chld1>yu)
                    {
                        chld1-=(yu-yl);
                    }
                    while(chld2<yl)
                    {
                        chld2+=(yu-yl);
                    }
                    while(chld2>yu)
                    {
                        chld2-=(yu-yl);
                    }*/
                    if(chld1<yl) chld1=rand()/(1.0*RAND_MAX);
                    if(chld2>yu) chld2=rand()/(1.0*RAND_MAX);
                    if (chld2 < yl) chld2 = rand()/(1.0*RAND_MAX);
                    if (chld2 > yu) chld2 = rand()/(1.0*RAND_MAX);
                }
                else
                {
                    chld1 = par1;
                    chld2 = par2;
                }
                newpop_ptr->ind[n].xreal[j] = chld1;
                newpop_ptr->ind[n+1].xreal[j] = chld2;
            }
        }
        else
        {
            for(j = 0;j < nvar;j++)
            {
                par1 = matepop_ptr->ind[y].xreal[j];
                par2 = matepop_ptr->ind[y+1].xreal[j]; 
                chld1 = par1;
                chld2 = par2;
                newpop_ptr->ind[n].xreal[j] = chld1;
                newpop_ptr->ind[n+1].xreal[j] = chld2;
            }
        }
        n = n+2; y=y+2;
    }
    return;
}

void bige_mutation(Population* newpop_ptr)
{
    int i,j,r;
    float rnd,delta,indi,*ptr,*ptr1,deltaq;
    float y,yl,yu,val,xy; 
    for(j = 0;j < popsize;j++)
    {
        for (i = 0;i < nvar; i++)
	    {
            rnd = rand()/(1.0*RAND_MAX);
	        //rnd=randomperc();
	        if(rnd <= pmut_r)
	        {
                y = newpop_ptr->ind[j].xreal[i];
	            yl = lim_r[i][0];
	            yu = lim_r[i][1]; 
	            if(y > yl)
       	        { 
                    if((y-yl) < (yu-y))
                        delta = (y - yl)/(yu - yl);
 	                else
  	                    delta = (yu - y)/(yu-yl);
                    rnd = rand()/(1.0*RAND_MAX); 
                    //rnd=randomperc();
                    indi = 1.0/(dim +1.0); 
                    if(rnd <= 0.5)
                    {
                        xy = 1.0-delta;
                        val = 2*rnd+(1-2*rnd)*(pow(xy,(dim+1)));
                        deltaq =  pow(val,indi) - 1.0;
                    }
                    else
                    {
                        xy = 1.0-delta;
                        val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(dim+1)));
                        deltaq = 1.0 - (pow(val,indi));
                    }
                    // Added by Deb (31/10/01)
                    y = y + deltaq * (yu-yl);
                    if (y < yl) y=yl; 
                    if (y > yu) y=yu;
                    newpop_ptr->ind[j].xreal[i] = y;
                }
                else // y == yl 
                {
                    xy = rand()/(1.0*RAND_MAX);
                    //xy=randomperc();
                    newpop_ptr->ind[j].xreal[i] = xy*(yu - yl) + yl;
                }
            }
            //  ptr++;
        }
    }
    return;
}

void bige_union_pop(Population* pop1,Population* pop2)
{
    int k=popsize;
    int i;
    for(k=popsize;k<2*popsize;k++)
    {
        for(i=0;i<nvar;i++)
        {
            (pop1->ind[k]).xreal[i]=(pop2->ind[k-popsize]).xreal[i];
        }
    }
}

void bige_nondominate_sort(Population* pop,LayerList* list)
{
    int nlayer=0;
    int i,j;
    Individual* indi,*indj;
    LayerIndSet* layer_ptr;
    DominIndSet* dominInd_ptr;
   // printf("yes\n");

    for(i=0;i<gpopsize;i++)
    {
        indi=&(pop->ind[i]);
        indi->dominated=0;
        for(j=0;j<gpopsize;j++)
        {
            if(j==i)
            {
                continue;
            }
            indj=&(pop->ind[j]);
            int comp=bige_comp_ind(indi,indj);
            if(comp>0)
            {
                bige_add_dominind(list,i,j);
               // indj->dominated++;
            }
            else if(comp<0)
            {
                indi->dominated++;
               // addDominInd(list,j,i);
            }
        }
        if(indi->dominated==0)
        {
            bige_add_layerind(list,i,nlayer);
            indi->rank=nlayer+1;
        }
    }
  //  printf("yes\n");
    layer_ptr=list->layers[nlayer];
    while(layer_ptr!=NULL)
    {
        nlayer++;
        LayerIndSet* temp=layer_ptr;
        while(temp!=NULL)
        {
            dominInd_ptr=list->dominset[temp->indNo];
            while(dominInd_ptr!=NULL)
            {
                int jj=dominInd_ptr->indNo;
                pop->ind[jj].dominated-=1;
                if(pop->ind[jj].dominated==0)
                {
                    pop->ind[jj].rank=nlayer+1;
                    bige_add_layerind(list,jj,nlayer);
                }
                dominInd_ptr=dominInd_ptr->next;
            }
            //printf("yes\n");
            temp=temp->next;
        }
        layer_ptr=list->layers[nlayer];
       // printf("yes\n");
    }
    list->nlayer=nlayer;
    return;
}

void bige_keepalive(Population* pop,LayerList* list, Population* nextpop_ptr)
{
    int i=0,j;
    int k=0;
    int indno;
    int inum=0;
    inum+=list->layerIndNum[i];
    LayerIndSet* layerind;
    Individual* indo;
    Individual* indm;
    while(inum<popsize)
    {
        layerind=list->layers[i];
        while(layerind!=NULL)
        {
            indno=layerind->indNo;
            indo=&(pop->ind[indno]);
            indm=&(nextpop_ptr->ind[k]);
            bige_copy_ind(indm,indo);
            layerind=layerind->next;
            k++;
        }
        i++;
        inum+=list->layerIndNum[i];
    }
    
    if(inum==popsize)
    {
        layerind=list->layers[i];
        while(layerind!=NULL)
        {
            indno=layerind->indNo;
            indo=&(pop->ind[indno]);
            indm=&(nextpop_ptr->ind[k]);
            bige_copy_ind(indm,indo);
            layerind=layerind->next;
            k++;
        }
    }
    else
    {
        int num=list->layerIndNum[i]; 
        inum-=num;
        layerind=list->layers[i];
        while(inum<=popsize)
        {
           /* int index=rand()%num;
            if(flags[index]!=0)
            {
                continue;
            }
            flags[index]=1;
            j=0;
            layerind=list->layers[i];
            while(j<index)
            {
                layerind=layerind->next;
                j++;
            }*/
            indno=layerind->indNo;
            indo=&(pop->ind[indno]);
            indm=&(nextpop_ptr->ind[k]);
            bige_copy_ind(indm,indo);
            k++;
            inum++;
            layerind=layerind->next;
        }
    }
    nextpop_ptr->maxrank=i+1;
}

void bige_clear_list(LayerList* list)
{
    LayerIndSet* temp1;
    DominIndSet* temp2;
    int i,j;
    for(i=0;i<gpopsize;i++)
    {
        temp2=list->dominset[i];
        while(temp2!=NULL)
        {
            DominIndSet* temp=temp2;
            temp2=temp2->next;
            free(temp);
        }
        list->dominset[i]=NULL;
    }
    for(j=0;j<list->nlayer;j++)
    {
        temp1=list->layers[j];
        while(temp1!=NULL)
        {
            LayerIndSet* temp=temp1;
            temp1=temp1->next;
            free(temp);
        }
        list->layers[j]=NULL;
    }
}

void bige_init_list(LayerList* list)
{
    bige_clear_list(list);
    memset(list->layerIndNum,0,sizeof(list->layerIndNum));
    list->nlayer=0;
    return;
}

void bige_add_layerind(LayerList* list,int i,int n)
{
    LayerIndSet* node=(LayerIndSet*)malloc(sizeof(LayerIndSet));
    node->indNo=i;
    node->next=list->layers[n];
    list->layers[n]=node;
    list->layerIndNum[n]+=1;
    return;
}

void bige_add_dominind(LayerList* list, int i,int j)
{
    DominIndSet* node=(DominIndSet*)malloc(sizeof(DominIndSet));
    node->indNo=j;
    node->next=list->dominset[i];
    list->dominset[i]=node;
    return;
}

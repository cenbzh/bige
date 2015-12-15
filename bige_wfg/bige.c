/*************************************************************************
	> File Name: bige.c
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年11月04日 星期三 19时38分42秒
 ************************************************************************/

#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bige.h"
#include <string.h>
#include <math.h>
#include <time.h>
void testout(Population*);
//extern float randomperc();
void bige_Engineer(FILE* outfile,char* problem,int type)
{
    Population* oldPop_ptr=&oldPops;
    Population* newPop_ptr=&newPops;
    Population* matePop_ptr=&matePops;
    Population* tempPop_ptr=&tempPops;
    LayerList* layerlist_ptr=&layerlist;
    printf("initPop..\n");
    initPop(oldPop_ptr);
    printf("initPop end\n");
    int i;
    radiu=pow(1.0/popsize,1.0/nfunc);
    for(i=0;i<ngener;i++)
    {
        printf("%d\n",i+1);
        //radiu=pow(1.0/popsize,1.0/nfunc);
        getObjectiveValue(oldPop_ptr,popsize,problem);
        proximityEstimation(oldPop_ptr,popsize);
        crowdingDegreeEstimation(oldPop_ptr,popsize);
        matingSelection(oldPop_ptr,matePop_ptr);
        generatingOffspring(matePop_ptr,newPop_ptr);
        getObjectiveValue(newPop_ptr,popsize,problem);
        initLayerList(layerlist_ptr);
        //radiu=pow(1.0/2*popsize,1.0/nfunc);
        environmentSelection(oldPop_ptr,newPop_ptr,matePop_ptr,layerlist_ptr,problem);
        copyPop(oldPop_ptr,matePop_ptr);
        if(i==ngener-1)
        {
            printf("\n\nRank Num:\n");
           // radiu=pow(1.0/popsize,1.0/nfunc);
           // outputRankNum(layerlist_ptr);
            printLayer(layerlist_ptr);
        }
        clearLayerList(layerlist_ptr);
    }
    getObjectiveValue(oldPop_ptr,popsize,problem);
    output(oldPop_ptr,outfile,type);
    for(i=0;i<nfunc;i++)
    {
        if(maxObjectives[i]<oldPop_ptr->maxObj[i])
        {
            maxObjectives[i]=oldPop_ptr->maxObj[i];
        }
    }
    return;
}

void initPop(Population* pop)
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
            d=rand()/(1.0*RAND_MAX);
            //d=randomperc();
            ind->xreal[j]=d;
        }
        //printf("\n");
        ind++;
    } 
    return;
}


void proximityEstimation(Population* pop,int size)
{
    int p,k;
    double sum=0;
    Individual* ind=&(pop->ind[0]);
    for(p=0;p<size;p++)
    {
        ind=&(pop->ind[p]);
        sum=0;
        for(k=0;k<nfunc;k++)
        {
            sum+=getNomalization(pop,p,k);
        }
        ind->proximity=sum;
    }
    return;
}

void crowdingDegreeEstimation(Population* pop,int size)
{
    getShareMatrix(pop,size);
    int p,k,q;
    Individual* ind=&(pop->ind[0]);
    for(p=0;p<size;p++)
    {
        ind->cr:owdingDegree=getCrowdingDegree(pop,p,size);
        ind++;
    }
    return;
}

double getCrowdingDegree(Population* pop,int p,int size)
{
    int q;
    double sum=0;
    for(q=0;q<size;q++)
    {
        if(q==p)
        {
            continue;
        }
        //double sh=shareFunction(pop,p,q);
        double sh=shMatrix[p][q];
        sum+=sh;
    }
    return pow(sum,1.0/2);
}

double getShareMatrix(Population* pop,int size)
{
    memset(shMatrix,0,sizeof(shMatrix));
    int i,j;
    Individual* indi,*indj;
    double dis=0;
    double Pri,Prj;
    double rnk;
    //printf("\n");
    for(i=0;i<size;i++)
    {
        indi=&(pop->ind[i]);
        for(j=i+1;j<size;j++)
        {
            indj=&(pop->ind[j]);
            dis=getDistance(pop,i,j);
            //printf("dis: %lf--- %lf\n",dis,radiu);
            if(dis>=radiu)
            {
                shMatrix[i][j]=0;
                shMatrix[j][i]=0;
            }
            else
            {
                //printf("Yes involve the shMatrix\n");
                Pri=indi->proximity;
                Prj=indj->proximity;
                if(Pri==Prj)
                {
                    rnk=rand()/(1.0*RAND_MAX);
                    if(rnk<=0.5)
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
                else if(Pri<Prj)
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
    /*for(i=0;i<popsize;i++)
    {
        for(j=0;j<popsize;j++)
        {
            printf("%12lf",shMatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");*/
    return 0;
}
double getDistance(Population* pop ,int p,int q)
{
    int k;
    double dis=0;
    for(k=0;k<nfunc;k++)
    {
        double fp=getNomalization(pop,p,k);
        double fq=getNomalization(pop,q,k);
        dis+=(fp-fq)*(fp-fq);
    }
    return pow(dis,1.0/2);
}

double getNomalization(Population* pop,int p,int k)
{
    if(pop->maxObj[k]==pop->minObj[k])
    {
        return 0;
    }
    Individual* ind=&(pop->ind[p]);
    return (ind->objs[k]-pop->minObj[k])/(pop->maxObj[k]-pop->minObj[k]);
}


void matingSelection(Population* oldPop,Population* matePop)
{
    int ind1,ind2;
    int p,k;
    Individual* select_ptr;
    Individual* ind1_ptr;
    Individual* ind2_ptr;
    Individual* mateInd_ptr;
    Individual* lastselect=NULL;
    for(p=0,k=0;p<popsize;)
    {
        totalcompInd++;
        ind1=rand()%popsize;
        ind2=rand()%popsize;
        while(ind2==ind1)
        {
            ind2=rand()%popsize;
        }
        ind1_ptr=&(oldPop->ind[ind1]);
        ind2_ptr=&(oldPop->ind[ind2]);
        mateInd_ptr=&(matePop->ind[k]);
        select_ptr=tournamentSelect(ind1_ptr,ind2_ptr);
        if(select_ptr==lastselect&&p%2==1)
        {
            nequal++;
            continue;
        }
        copyInd(mateInd_ptr,select_ptr);
        lastselect=select_ptr;
        p++;
        k++;
    }
    return;
}
/*
 * 在NsgaII 中它的环境选择利用了rank来比较的，这样就会考虑到了个体在种群位置的信息，如果直接单独地比较两个个体的话，就没有考虑到个体在种群位置的其他信息
* */
Individual* tournamentSelect(Individual* ind1_ptr,Individual* ind2_ptr)
{
    int r=compInd(ind1_ptr,ind2_ptr);
    if(r>0)
    {
        return ind1_ptr;
    }
    else if(r<0)
    {
        return ind2_ptr;
    }
    else
    {
        double pr=rand()/(1.0*RAND_MAX);
        if(pr<=0.5)
        {
            return ind1_ptr;
        }
        else
        {
            return ind2_ptr;
        }
    }
}

int compInd(Individual* ind1_ptr,Individual* ind2_ptr)
{
    //totalcompInd++;
    double pr1,pr2,cd1,cd2;
    pr1=ind1_ptr->proximity;
    pr2=ind2_ptr->proximity;
    cd1=ind1_ptr->crowdingDegree;
    cd2=ind2_ptr->crowdingDegree;
    /*if(fabs(pr1-pr2)>0.000001)
    {
        if(pr1<pr2)
        {
            if(fabs(cd1-cd2)>0.000001&&cd1<cd2)
            {
                return 1;
            }
        }
        else if(pr2<pr1)
        {
            if(fabs(cd1-cd2)>0.000001&&cd2<cd1)
            {
                return -1;
            }
        }
        nequal++;
        return 0;
    }
    else
    {
        nequal++;
        return 0;
    }*/
    if((pr1<pr2&&cd1<cd2)||(pr1<pr2&&cd1<cd2))
    {
        return 1;
    }
    if((pr1>pr2&&cd1>cd2))
    {
        return -1;
    }
    if(pr1==pr2||cd1==cd2)
    {
        //nequal++;
    }
    
    return 0;
}
void copyInd(Individual* desInd,Individual* scrInd)
{
    desInd->rank=scrInd->rank;
    desInd->dominated=scrInd->dominated;
    int i;
    for(i=0;i<nvar;i++)
    {
        desInd->xreal[i]=scrInd->xreal[i];
    }
    for(i=0;i<nfunc;i++)
    {
        desInd->objs[i]=scrInd->objs[i];
    }
    desInd->proximity=scrInd->proximity;
    desInd->crowdingDegree=scrInd->crowdingDegree;
    return ;
}


void generatingOffspring(Population* matePop,Population* newPop)
{
    crossover(matePop,newPop);
    mutation(newPop);
    return ;
}

void crossover(Population* matePop,Population* newPop)
{
    //extern float randomperc();
   /* int i,j,y,n,r;
    double rnd,par1,par2,chld1,chld2,betaq,beta,alpha;
    double y1,y2,yu,yl,expp;
    y=0;
    n=0;
    for(i=0;i<popsize/2;i++)
    {
        rnd=rand()/(1.0*RAND_MAX);
        if(rnd<=pcross)
        {
            for(j=0;j<nvar;j++)
            {
                par1=matePop->ind[y].xreal[j];
                par2=matePop->ind[y+1].xreal[j];
                yl=lim_r[j][0];
                yu=lim_r[j][1];
                rnd=rand()/(1.0*RAND_MAX);
                if(rnd<=0.5)
                {
                    if(fabs(par1-par2)>0.000001)
                    {
                        if(par2>par1)
                        {
                            y2=par2;
                            y1=par1;
                        }
                        else{
                            y1=par2;
                            y2=par1;
                        }
                        if((y1-yl)>(yu-y2))
                        {
                            beta=1+(2*(yu-y2)/(y2-y1));
                        }
                        else
                        {
                            beta=1+(2*(y1-yl)/(y2-y1));
                        }
                        expp=di+1.0;
                        beta=1.0/beta;
                        alpha=2.0-pow(beta,expp);
                        if(alpha<0.0)
                        {
                            printf("ERROR!!!\n");
                            exit(-1);
                        }
                        rnd=rand()/(1.0*RAND_MAX);
                        if(rnd<=1.0/alpha)
                        {
                            alpha=alpha*rnd;
                            expp=1.0/(di+1.0);
                            betaq=pow(alpha,expp);
                        }
                        else
                        {
                            alpha=alpha*rnd;
                            alpha=1.0/(2.0-alpha);
                            expp=1.0/(di+1.0);
                            if(alpha<0.0)
                            {
                                printf("ERROR!!!!\n");
                                exit(-1);
                            }
                            betaq=pow(alpha,expp);
                        }
                        chld1=0.5*(y1+y2)-betaq*(y2-y1);
                        chld2=0.5*(y1+y2)+betaq*(y2-y1);
                    }
                    else{
                        betaq=1.0;
                        y1=par1;y2=par2;
                        chld1=0.5*(y1+y2)-betaq*(y2-y1);
                        chld2=0.5*(y1+y2)+betaq*(y2-y1);
                    }
                    if(chld1<yl) chld1=yl;
                    if(chld2>yu) chld2=yu;
                    if(chld2<yl) chld2=yl;
                    if(chld1>yu) chld1=yu;
                }
                else
                {
                    chld1=par1;
                    chld2=par2;
                }
                newPop->ind[n].xreal[j]=chld1;
                newPop->ind[n+1].xreal[j]=chld2;
            }
        }
        else
        {
            for(j=0;j<nvar;j++)
            {
                par1=matePop->ind[y].xreal[j];
                par2=matePop->ind[y+1].xreal[j];
                chld1=par1;
                chld2=par2;
                newPop->ind[n].xreal[j]=chld1;
                newPop->ind[n+1].xreal[j]=chld2;
            }
        }
        n=n+2;
        y=y+2;
    }
    return;*/
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
            par1 = matePop->ind[y].xreal[j];
            par2 = matePop->ind[y+1].xreal[j]; 
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
                newPop->ind[n].xreal[j] = chld1;
                newPop->ind[n+1].xreal[j] = chld2;
            }
        }
        else
        {
            for(j = 0;j < nvar;j++)
            {
                par1 = matePop->ind[y].xreal[j];
                par2 = matePop->ind[y+1].xreal[j]; 
                chld1 = par1;
                chld2 = par2;
                newPop->ind[n].xreal[j] = chld1;
                newPop->ind[n+1].xreal[j] = chld2;
            }
        }
        n = n+2; y=y+2;
    }
    return;
}


void mutation(Population* newPop)
{
    //extern float randomperc();
    /*int i,j,r;
    double rnd,delta,indi,*ptr,*ptr1,deltaq;
    double y,yl,yu,val,xy;
    for(j=0;j<popsize;j++)
    {
        for(i=0;i<nvar;i++)
        {
            rnd=rand()/(1.0*RAND_MAX);
            if(rnd<=pmut_r)
            {
                y=newPop->ind[j].xreal[i];
                yl=lim_r[i][0];
                yu=lim_r[i][1];
                if(y>yl)
                {
                    if((y-yl)<(yu-y))
                        delta=(y-yl)/(yu-yl);
                    else
                        delta=(yu-y)/(yu-yl);
                    rnd=rand()/(1.0*RAND_MAX);
                    indi=1.0/(dim+1.0);
                    if(rnd<=0.5)
                    {
                        xy=1.0-delta;
                        val=2*rnd+(1-2*rnd)*(pow(xy,dim+1));
                        deltaq=pow(val,indi)-1.0;
                    }
                    else
                    {
                        xy=1.0-delta;
                        val=2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,dim+1));
                        deltaq=1.0-pow(val,indi);
                    }
                    y=y+deltaq*(yu-yl);
                    if(y<yl) y=yl;
                    if(y>yu) y=yu;
                    newPop->ind[j].xreal[i]=y;
                }
                else
                {
                    xy=rand()/(1.0*RAND_MAX);
                    newPop->ind[j].xreal[i]=xy*(yu-yl)+yl;
                }
            }
        }
    }
    return;*/
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
	      y = newPop->ind[j].xreal[i];
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
		  newPop->ind[j].xreal[i] = y;
		}
	      else // y == yl 
		{
		  xy = rand()/(1.0*RAND_MAX);
		  //xy=randomperc();
            newPop->ind[j].xreal[i] = xy*(yu - yl) + yl;
		}
	    }
	  //  ptr++;
	}
    }
  return ;
    
}

void unionPop(Population* pop1,Population* pop2)
{
    int k=popsize;
    int i;
    for(k=popsize;k<2*popsize;k++)
    {
        for(i=0;i<nvar;i++)
        (pop1->ind[k]).xreal[i]=(pop2->ind[k-popsize]).xreal[i];
    }
}

void environmentSelection(Population* oldPop,Population* newPop,Population*matePop,LayerList* layerlist_ptr,char* problem)
{
    unionPop(oldPop,newPop);
    getObjectiveValue(oldPop,2*popsize,problem);
    /*gpopsize=deleteDuplicate(oldPop,2*popsize);
    if(gpopsize<popsize)
    {
        return;
    }*/
    gpopsize=2*popsize;
    proximityEstimation(oldPop,gpopsize);
    crowdingDegreeEstimation(oldPop,gpopsize);
    //printf("yes\n");
    nondominatedSort(oldPop,layerlist_ptr);
   // printf("yes\n");
    keepalive(oldPop,layerlist_ptr,matePop);
    return ;
}

void nondominatedSort(Population* oldPop,LayerList* list)
{
    int nlayer=0;
    int i,j;
    Individual* indi,*indj;
    LayerInd* layer_ptr;
    DominInd* dominInd_ptr;
   // printf("yes\n");

    for(i=0;i<gpopsize;i++)
    {
        indi=&(oldPop->ind[i]);
        indi->dominated=0;
        for(j=0;j<gpopsize;j++)
        {
            if(j==i)
            {
                continue;
            }
            indj=&(oldPop->ind[j]);
            int comp=compInd(indi,indj);
            if(comp>0)
            {
                addDominInd(list,i,j);
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
            addLayerInd(list,i,nlayer);
            indi->rank=nlayer+1;
        }
    }
  //  printf("yes\n");
    layer_ptr=list->layers[nlayer];
    while(layer_ptr!=NULL)
    {
        nlayer++;
        LayerInd* temp=layer_ptr;
        while(temp!=NULL)
        {
            dominInd_ptr=list->dominlist[temp->i];
            while(dominInd_ptr!=NULL)
            {
                int jj=dominInd_ptr->j;
                oldPop->ind[jj].dominated-=1;
                if(oldPop->ind[jj].dominated==0)
                {
                    oldPop->ind[jj].rank=nlayer+1;
                    addLayerInd(list,jj,nlayer);
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

void addLayerInd(LayerList* list,int i,int layer_no)//i->layer
{
    LayerInd* node=(LayerInd*)malloc(sizeof(LayerInd));
    node->i=i;
    node->next=list->layers[layer_no];
    list->layers[layer_no]=node;
    list->layerNum[layer_no]+=1;
    return;
}
void addDominInd(LayerList* list,int i,int j)//j->i i dominate j
{
    DominInd* node=(DominInd*)malloc(sizeof(DominInd));
    node->j=j;
    node->next=list->dominlist[i];
    list->dominlist[i]=node;
    return;
}

void initLayerList(LayerList* list)
{
    clearLayerList(list);
    memset(list->layerNum,0,sizeof(list->layerNum));
    list->nlayer=0;
    return;
}

void clearLayerList(LayerList* list)
{
    LayerInd* temp1;
    DominInd* temp2;
    int i,j;
    for(i=0;i<gpopsize;i++)
    {
        temp2=list->dominlist[i];
        while(temp2!=NULL)
        {
            DominInd* temp=temp2;
            temp2=temp2->next;
            free(temp);
        }
        list->dominlist[i]=NULL;
    }
    for(j=0;j<list->nlayer;j++)
    {
        temp1=list->layers[j];
        while(temp1!=NULL)
        {
            LayerInd* temp=temp1;
            temp1=temp1->next;
            free(temp);
        }
        list->layers[j]=NULL;
    }
    
}

void keepalive(Population* oldPop,LayerList* list,Population* matePop)
{
    int i=0,j;
    int k=0;
    int indno;
    int inum=0;
    inum+=list->layerNum[i];
    LayerInd* layerind;
    Individual* indo;
    Individual* indm;
    while(inum<popsize)
    {
        layerind=list->layers[i];
        while(layerind!=NULL)
        {
            indno=layerind->i;
            indo=&(oldPop->ind[indno]);
            indm=&(matePop->ind[k]);
            copyInd(indm,indo);
            layerind=layerind->next;
            k++;
        }
        i++;
        inum+=list->layerNum[i];
    }
    
    if(inum==popsize)
    {
        layerind=list->layers[i];
        while(layerind!=NULL)
        {
            indno=layerind->i;
            indo=&(oldPop->ind[indno]);
            indm=&(matePop->ind[k]);
            copyInd(indm,indo);
            layerind=layerind->next;
            k++;
        }
    }
    else
    {
        int num=list->layerNum[i];
        
  
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
            indno=layerind->i;
            indo=&(oldPop->ind[indno]);
            indm=&(matePop->ind[k]);
            copyInd(indm,indo);
            k++;
            inum++;
            layerind=layerind->next;
        }
    }
    matePop->maxrank=i+1;
}

void copyPop(Population* desPop,Population* scrPop)
{
    int p;
    Individual* ind1;
    Individual* ind2;
    desPop->maxrank=scrPop->maxrank;
    for(p=0;p<popsize;p++)
    {
        ind1=&(desPop->ind[p]);
        ind2=&(scrPop->ind[p]);
        copyInd(ind1,ind2);
    }
    int i;
    for(i=0;i<nfunc;i++)
    {
        desPop->maxObj[i]=scrPop->maxObj[i];
        desPop->minObj[i]=scrPop->minObj[i];
    }
    return ;
}


void getObjectiveValue(Population* pop,int size,char* problem)
{
    int p;
    int j;
    Individual* ind;
    //printf("%d\n",size);
    for(p=0;p<size;p++)
    {
        //printf("yes %d\n",p);
        ind=&(pop->ind[p]);
        wfg_eval(ind->xreal,nvar,2*(nfunc-1),nfunc,problem,ind->objs);
        //DTLZ(ind->xreal,ind->objs,1,nvar,nfunc);
        for(j=0;j<nfunc;j++)
        {
            if(p==0)
            {
                pop->maxObj[j]=ind->objs[j];
                pop->minObj[j]=ind->objs[j];
            }
            else
            {
                if(pop->maxObj[j]<ind->objs[j])
                {
                    pop->maxObj[j]=ind->objs[j];
                }
                if(pop->minObj[j]>=ind->objs[j])
                {
                    pop->minObj[j]=ind->objs[j];
                }
            }
            //pop->maxObj[j]=0.5;
            //pop->minObj[j]=0;
        }
        /*double sum=0;
        for(j=0;j<nfunc;++j)
        {
        	sum+=ind->objs[j];
        }
        	if(min>sum)
        	min=sum;*/
    }
    /*for(j=0;j<nfunc;++j)
    {
    pop->maxObj[j]=min/nfunc;
    pop->minObj[j]=0;
    }*/
}

void input(FILE* infile)
{
    int generNum;
    int varNum;
    int funcNum;
    double crossP;
    double mutP;
    double lim_b;
    double lim_u;
    int popNum;
    double ditemp,dimtemp;
    int ans;
    int i;
    printf("Now Input the parameters of BIGE:\n");
    printf("No of generations: ");
    fscanf(infile,"%d",&generNum);
    printf("The size of Population(<=%d): ",MAXPOP);
    fscanf(infile,"%d",&popNum);
    printf("The number of variables(<=%d): ",MAXVAR);
    fscanf(infile,"%d",&varNum);
    nvar=varNum;
    printf("The number of objectives(<=%d): ",MAXFUN);
    fscanf(infile,"%d",&funcNum);
    printf("The probability of crossover: ");
    fscanf(infile,"%lf",&crossP);
    printf("The probability of mutation: ");
    fscanf(infile,"%lf",&mutP);
    printf("Are there scope for variables(0 or 1):");
    fscanf(infile,"%d",&ans);
    printf("%d\n",ans);
    if(ans==1)
    {
       // printf("Input the scope of variables:\n");
        for(i=0;i<nvar;i++)
        {
            //fscanf(infile,"%lf%lf",&lim_b,&lim_u);
            //lim_r[i][0]=lim_b;
            //lim_r[i][1]=lim_u;
            lim_r[i][0]=0;
            lim_r[i][1]=1;
            //printf("%lf %lf\n",lim_r[i][0],lim_r[i][1]);
        }
    }
    printf("Input the Distribution index of crossover: ");
    fscanf(infile,"%lf",&ditemp);
    printf("Input the Distribution index of mutation: ");
    fscanf(infile,"%lf",&dimtemp);
    ngener=generNum;
    popsize=popNum;
    //nvar=varNum;
    nfunc=funcNum;
    pcross=crossP;
    //pmut_r=mutP;
    pmut_r=1.0/(nvar);
    isLimit=ans;
    di=ditemp;
    dim=dimtemp;
    printf("\n\n");
    return;
}

void output(Population* pop,FILE* outfile,int type)
{
    //printf("The maxrank is %d\n",pop->maxrank);
    int i,j;
    Individual* ind;
    //printf("%6s%12s%12s%12s%12s\n","No.","Pr.","Cd.","Rank.","objs.");
    int count=popsize;
    /*if(type==1)
    {
    
    }
    else
    {
        int k=0;
        memset(dupFlags,0,sizeof(dupFlags));
        Individual* indi;
        Individual* indj;
        for(i=0;i<popsize;i++)
        {
            indi=&(pop->ind[i]);
            for(j=0;j<popsize;j++)
            {
                if(j==i)
                {
                    continue;
                }
                indj=&(pop->ind[j]);
                int p=0;
                for(p=0;p<nvar;p++)
                {
                    if(indi->xreal[p]<indj->xreal[p])
                    {
                        break;
                    }
                }
                if(p==nvar)
                {
                    break;
                }
            }
            if(j==popsize)
            {
                k++;
                dupFlags[i]=1;
            }
        }
        count=k;
        
    }*/
    //fprintf(outfile,"%12d\n",count);
    int printflag=1;
    for(i=0;i<popsize;i++)
    {
       /* if(type==0&&dupFlags[i]==0)
        {
           printflag=0;
        }
        else
        {
           printflag=1;
        }*/
        if(printflag==1)
        {
        ind=&(pop->ind[i]);
        //printf("%6d%12.4lf%12.4lf%12d",i+1,ind->proximity,ind->crowdingDegree,ind->rank);
        for(j=0;j<nfunc;j++)
        {
            fprintf(outfile,"%12lf",ind->objs[j]);
        }
        fprintf(outfile,"\n");
        }
       /* for(j=0;j<nvar;j++)
        {
            printf("%12.4lf",ind->xreal[j]);
        }
        printf("\n---------------------------\n");*/
    }
    return ;
}

int main(int argc,char* argv[])
{
    if(argc!=2)
    {
        printf("usage: %s problem\n",argv[0]);
        exit(-1);
    }
    extern void warmup_random(float);
    int i;
    srand(time(NULL));
    FILE* infile;
    FILE* outfile;
    char inName[200];
    strcpy(inName,DIR_PATH);
    strcat(inName,argv[1]);
    strcat(inName,"/");
    strcat(inName,argv[1]);
    strcat(inName,"_input.txt");
    printf("%s\n",inName);
    //infile=fopen(inName,"r");
    char outName[200];
    char format[10];
    int type=1;
    if(strcmp(argv[1],"WFG1")==0||strcmp(argv[1],"WFG2")==0||strcmp(argv[1],"WFG3")==0)
    {
        type=0;
    }
    memset(maxObjectives,0,sizeof(maxObjectives));
    int runtime=1;
    warmup_random(0.4);
    for(i=0;i<runtime;i++)
    {
        nequal=0;
        totalcompInd=0;
        memset(outName,'\0',sizeof(outName));
        memset(format,'\0',sizeof(format));
        infile=fopen(inName,"r");
        input(infile);
        strcpy(outName,DIR_PATH);
        strcat(outName,argv[1]);
        strcat(outName,"/finalFit2_");
        itoa(i+1,format,10);
        strcat(outName,format);
        outfile=fopen(outName,"w");
        bige_Engineer(outfile,argv[1],type);
       // fclose(infile);
        if(i==runtime-1)
        {
            int j;
            for(j=0;j<nfunc;j++)
            {
                fprintf(outfile,"%12lf",maxObjectives[j]);
            }
        }
        fclose(outfile);
        fclose(infile);
        printf("T:%d  E:%d\n",totalcompInd,nequal);
    }
    //fclose(infile);
    return 0;
    //output()
}

/*
void testout(Population* pop)
{
    int p;
    int v;
    Individual* ind;
    for(p=0;p<popsize;p++)
    {
        ind=&(pop->ind[p]);
        for(v=0;v<nvar;v++)
        {
            printf("%lf ",ind->xreal[v]);
        }
        printf("\n");
    }
    printf("\n");
}

void outputRankNum(LayerList* list)
{
    int num=list->nlayer;
    int i;
    printf("%12s%12s%12s:%d\n","LayerNo.","Num.","R.",gpopsize);
    for(i=0;i<num;i++)
    {
        printf("%12d%12d%11.2f%s\n",i+1,list->layerNum[i],list->layerNum[i]/(1.0*gpopsize)*100,"%");
    }
    return ;
}

double dominatedRate(Population* finalPop,Population* initPop,double *rates)
{
    int sum1=0;
    int sum2=0;
    int flag=0;
    int i,j;
    Individual* find;
    Individual* iind;
    for(i=0;i<popsize;i++)
    {
        find=&(finalPop->ind[i]);
        iind=&(initPop->ind[i]);
        for(j=0;j<nfunc;j++)
        {
            if(find->objs[j]>iind->objs[j])
            {
                break;
            }
        }
        if(j==nfunc)
        {
            sum1++;
        }
        for(j=0;j<nfunc;j++)
        {
            if(find->objs[j]<iind->objs[j])
            {
                break;
            }
        }
        if(j==nfunc)
        {
            sum2++;
        }
    }
    rates[0]=sum1/(1.0*popsize);
    rates[1]=(popsize-sum1-sum2)/(1.0*popsize);
    rates[2]=sum2/(1.0*popsize);
    return 0;
}

int deleteDuplicate(Population* pop,int size)
{
    memset(dupFlags,0,sizeof(dupFlags));
    int i,j;
    Individual* indi;
    Individual* indj;
    for(i=0;i<size;i++)
    {
        indi=&(pop->ind[i]);
        if(dupFlags[i]==0)
        {
            for(j=i+1;j<size;j++)
            {
                indj=&(pop->ind[j]);
                if(isSameValue(indi,indj)==1)
                {
                    dupFlags[j]=1;
                }
            }
        }
    }
    Individual* last=&(pop->ind[0]);
    Individual* indc;
    int result=0;
    for(i=0;i<size;i++)
    {
        indc=&(pop->ind[i]);
        if(dupFlags[i]==0)
        {
            copyInd(last,indc);
            last++;
            result++;
        }
    }
    return result;
}

int isSameValue(Individual* ind1,Individual* ind2)
{
    int k;
    for(k=0;k<nfunc;k++)
    {
        if(ind1->objs[k]!=ind2->objs[k])
        {
            return 0;
        }
    }
    /*int j;
    for(j=0;j<nvar;j++)
    {
        if(ind1->xreal[j]!=ind2->xreal[j])
        {
            return 0;
        }
    }
    return 1;
}
*/
void printLayer(LayerList* list)
{
    int n=list->nlayer;
    int i;
    printf("\nLayer: \n");
    for(i=0;i<n;i++)
    {
        printf("Layer %d: ",i+1);
        LayerInd* layerInd=list->layers[i];
        while(layerInd!=NULL)
        {
            printf("%d ",layerInd->i);
            layerInd=layerInd->next;
        }
        printf("\n");
    }
    printf("========\n");
}

void itoa(int num,char str[],int n)
{
    int i=0;
    int sign=1;
    if(num<0)
    {
        sign=-1;
    }
    int temp=num;
    int count=0;
    while(temp!=0)
    {
        count++;
        temp/=10;
    }
    if(sign==-1)
    {
        str[i]='-';
        i++;
        count++;
    }
    if(count>n)
    {
        printf("error\n");
        exit(-1);
    }
    for(i=count-1;i>=0;i--)
    {
        str[i]=num%10+'0';
        num=num/10;
        if(num==0)
        {
            break;
        }
    }
    
}

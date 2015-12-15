/*************************************************************************
	> File Name: bige.h
	> Author: Binzhong Cen
	> Mail: cenbzh@gmail.com
	> Created Time: 2015年11月02日 星期一 17时58分37秒
 ************************************************************************/

#ifndef _BIGE_H
#define _BIGE_H
#define MAXPOP 1000
#define MAXVAR 500
#define MAXFUN 30
#define MAXCONS 20
#define DIR_PATH "/home/cenbzh/Documents/paper/nsgaii/c/obj5/"

int nvar;   /*决策变量个数*/
int ngener;/*遗传代数*/
int isLimit;/*决策变量是否有范围*/
double pcross;/*遗传交叉概率*/
double pmut_r;/*遗传突变概率*/
double lim_r[MAXVAR][2];/*决策变量的范围*/
int nfunc;/*目标个数*/
int popsize;/*种群大小*/
double radiu;/*小生境半径*/
double di,dim;/*交叉和遗传要用到的分布指标*/
typedef struct Individual/*解个体的结构*/
{
    int rank;/*记录个体的rank值*/
    int dominated;/*记录个体被其他个体个体占优的个数*/
    double xreal[MAXVAR];/*记录决策变量的值*/
    double objs[MAXFUN];/*记录目标函数的值*/
    double proximity;/*收敛度*/
    double crowdingDegree;/*拥挤度*/
}Individual;

typedef struct Population/*种群的结构*/
{
    int maxrank;/*分层排序后的最大rank值*/
    double maxObj[MAXFUN];/*种群在每个目标函数上的最大值*/
    double minObj[MAXFUN];/*种群在每个目标函数上的最小值*/
    Individual ind[2*MAXPOP];/*种群的个体*/
}Population;



typedef struct LayerInd/*用单链表来存放同一rank上的个体*/
{
    int i;/*个体在种群的下标*/
    struct LayerInd* next;
}LayerInd;

typedef struct DominInd/*用链表来存放被该个体占优的个体*/
{
    int j;/*被占优个体对应在种群的下标*/
    struct DominInd* next;
}DominInd;

typedef struct LayerList
{
    int nlayer;/*分层排序后得到的层数*/
    LayerInd* layers[2*MAXPOP];
    DominInd* dominlist[2*MAXPOP];
    int layerNum[2*MAXPOP];//从0开始的
}LayerList;
Population oldPops;
Population matePops;
Population newPops;
LayerList layerlist;
Population tempPops;
double shMatrix[2*MAXPOP][2*MAXPOP];/*记录共享函数的值*/
int dupFlags[2*MAXPOP];
int gpopsize;
double maxObjectives[MAXFUN];


int totalcompInd;
int nequal;

extern int wfg_eval(double*,int ,int ,int ,char*,double*);/*获得wfg问题的目标函数值*/
extern void DTLZ(double x[],double f[],int which,int dim,int nobj);

void bige_Engineer(FILE*,char*,int);/*bige的主框架*/

void proximityEstimation(Population* pop,int size);/*收敛度的估算*/

void crowdingDegreeEstimation(Population* pop,int size);/*拥挤度的估算*/

void getObjectiveValue(Population* pop,int size,char* problem);/*获得目标函数的值*/

void initPop(Population* pop);/*初始化种群*/

void matingSelection(Population* oldPop,Population* matePop);/*繁殖选择*/

void generatingOffspring(Population* matePop,Population* newPop);/*产生后代*/

void environmentSelection(Population* oldPop,Population* newPop,Population* nextPop,LayerList*,char* );/*环境选择*/

Individual* tournamentSelect(Individual* d1,Individual* d2);/*锦标赛选择*/

void copyPop(Population* desPop,Population* scrPop);/*复制种群*/

double getCrowdingDegree(Population* pop,int p,int size);/*获得个体p的拥挤度*/

//double shareFunction(Population* pop,int p,int q);/*共享函数的计算*/

double getDistance(Population* pop,int p,int q);/*欧式距离的计算*/

double getNomalization(Population* pop,int p,int k);/*归一化*/

int compInd(Individual*,Individual* );/*判断两个个体的占优关系*/

void copyInd(Individual* desInd,Individual* scrInd);/*复制个体*/

void crossover(Population*,Population*);/*交叉操作*/

void mutation(Population*);/*突变操作*/

void unionPop(Population*,Population*);/*合并两个种群*/

void nondominatedSort(Population*,LayerList* list);/*非占优排序*/

void keepalive(Population*,LayerList*,Population*);/*选出新的进化种群*/

void clearLayerList(LayerList*);/*清理LayerList的结构中的内容*/
void initLayerList(LayerList*);/*初始化LayerList的结构*/
void addLayerInd(LayerList*,int ,int );/*将个体添加到相应的分层*/
void addDominInd(LayerList*,int ,int );/*将个体添加到某个个体的占优链表中*/

void input(FILE*);/*输入函数*/
void output(Population* pop,FILE*,int);/*输出函数*/
/*void outputRankNum(LayerList* );

double dominatedRate(Population* finalPop,Population* initPop,double* rates);
*/
double getShareMatrix(Population*,int);
/*int deleteDuplicate(Population*,int);
int isSameValue(Individual*,Individual*);*/
void printLayer(LayerList*);
void itoa(int num,char str[],int n);
#endif

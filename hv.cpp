/*

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

FILECONTENTS *readFile(char[]);

extern void printContents(FILECONTENTS *);

//-----------------read------------------
static void trimLine(char line[])
{
	int i = 0;

	while(line[i] != '\0')
	{
		if (line[i] == '\r' || line[i] == '\n')
		{
			line[i] = '\0';
			break;
		}
		i++;
	}
}

void printContents(FILECONTENTS *f)
{
	for (int i = 0; i < f->nFronts; i++)
	{
		printf("Front %d:\n", i+1);
		for (int j = 0; j < f->fronts[i].nPoints; j++)
		{
			printf("\t");
			for (int k = 0; k < f->fronts[i].n; k++)
			{
				printf("%f ", f->fronts[i].points[j].objectives[k]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

int dim;
POINT ref; // the reference point 
FILECONTENTS *readFile(char filename[])
{ 	
	FILE *fp;
	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		printf("File could not be opened\n");
		exit(1);
	}
	FILECONTENTS *fc = (FILECONTENTS *)malloc(sizeof(FILECONTENTS));
	fc->nFronts = 30;		// number of pfs
	fc->fronts = (FRONT*)malloc(sizeof(FRONT) * fc->nFronts);
  	for (int i = 0; i < fc->nFronts; i++) 
    { 
    	//fscanf(fp,"%d",&fc->fronts[i].nPoints);
    	fc->fronts[i].nPoints=100;
    	fc->fronts[i].n = dim;
    	fc->fronts[i].points = (POINT*)malloc(sizeof(POINT) * fc->fronts[i].nPoints);
    	for(int j = 0; j < fc->fronts[i].nPoints; ++j)
    	{
    		fc->fronts[i].points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * dim);
    		for(int k = 0; k < dim; ++k)
    		{
    			fscanf(fp,"%lf",&fc->fronts[i].points[j].objectives[k]);
    			if(fc->fronts[i].points[j].objectives[k]>ref.objectives[k])
    			{
    				fc->fronts[i].points[j].objectives[k]=ref.objectives[k];
   				}
   			}
   		}
    }
    return fc;
}
/*
FILECONTENTS *readFile(char filename[])
{
	FILE *fp;
	char line[BUFSIZ];
	int front = 0, point = 0, objective = 0;

	FILECONTENTS *fc = (FILECONTENTS *)malloc(sizeof(FILECONTENTS));
	fc->nFronts = 0;
	fc->fronts = NULL;

	fp = fopen(filename, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "File %s could not be opened\n", filename);
		exit(EXIT_FAILURE);
	}

	while(fgets(line, sizeof line, fp) != NULL)
	{
		trimLine(line);
		if (strcmp(line, "#") == 0)
		{
			front = fc->nFronts;
			fc->nFronts++;
			fc->fronts = (FRONT *)realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
			fc->fronts[front].nPoints = 0;
			fc->fronts[front].points = NULL;
		}
		else
		{
			FRONT *f = &fc->fronts[front];
			point = f->nPoints;
			f->nPoints++;
			f->points = (POINT*)realloc(f->points, sizeof(POINT) * f->nPoints);
			f->n = 0;
			f->points[point].objectives = NULL;
			char *tok = strtok(line, " \t\n");
			do
			{
				POINT *p = &f->points[point];
				objective = f->n;
				f->n++;
				p->objectives = (OBJECTIVE*)realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
				p->objectives[objective] = atof(tok);
			} while ((tok = strtok(NULL, " \t\n")) != NULL);
		}
	}

	fc->nFronts--;
	// for (int i = 0; i < fc->nFronts; i++) fc->fronts[i].n = fc->fronts[i].points[0].nObjectives;
        fclose(fp);
	//printf("Read %d fronts\n", fc->nFronts);
	   //printContents(fc); 
	return fc;
}
//-------------------------read end------------------------------------------
*/


#define BEATS(x,y)   (x > y) 
#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y)) 

int n;     // the number of objectives 


FRONT *fs;    // memory management stuff 
int fr = 0;   // current depth 
int maxm = 0; // identify the biggest fronts in the file 
int maxn = 0; 
int safe;     // the number of points that don't need sorting 

double totaltime;

double hv(FRONT);


int greater(const void *v1, const void *v2)
// this sorts points worsening in the last objective
{
  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  for (int i = n - 1; i >= 0; i--)
    if BEATS(p.objectives[i],q.objectives[i]) return -1;
    else
    if BEATS(q.objectives[i],p.objectives[i]) return  1;
  return 0;
}


int greaterabbrev(const void *v1, const void *v2)
// this sorts points worsening in the penultimate objective
{
  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  for (int i = n - 2; i >= 0; i--)
    if BEATS(p.objectives[i],q.objectives[i]) return -1;
    else
    if BEATS(q.objectives[i],p.objectives[i]) return  1;
  return 0;
}


int dominates2way(POINT p, POINT q, int k)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
// k is the highest index inspected 
{
  for (int i = k; i >= 0; i--)
    if BEATS(p.objectives[i],q.objectives[i]) 
      {for (int j = i - 1; j >= 0; j--) 
         if BEATS(q.objectives[j],p.objectives[j]) return 0; 
       return -1;}
    else
    if BEATS(q.objectives[i],p.objectives[i]) 
      {for (int j = i - 1; j >= 0; j--) 
         if BEATS(p.objectives[j],q.objectives[j]) return 0; 
       return  1;}
  return 2;
}


bool dominates1way(POINT p, POINT q, int k)
// returns true if p dominates q or p == q, false otherwise 
// the assumption is that q doesn't dominate p 
// k is the highest index inspected 
{
  for (int i = k; i >= 0; i--)
    if BEATS(q.objectives[i],p.objectives[i]) 
      return false;
  return true;
}


void makeDominatedBit(FRONT ps, int p)
// creates the front ps[0 .. p-1] in fs[fr], with each point bounded by ps[p] and dominated points removed 
{
  int l = 0;
  int u = p - 1;
  for (int i = p - 1; i >= 0; i--)
    if (BEATS(ps.points[p].objectives[n - 1],ps.points[i].objectives[n - 1]))
      {fs[fr].points[u].objectives[n - 1] = ps.points[i].objectives[n - 1]; 
       for (int j = 0; j < n - 1; j++) 
	 fs[fr].points[u].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
       u--;}
    else
      {fs[fr].points[l].objectives[n - 1] = ps.points[p].objectives[n - 1]; 
       for (int j = 0; j < n - 1; j++) 
 	 fs[fr].points[l].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
       l++;}
  POINT t;
  // points below l are all equal in the last objective; points above l are all worse 
  // points below l can dominate each other, and we don't need to compare the last objective 
  // points above l cannot dominate points that start below l, and we don't need to compare the last objective 
  fs[fr].nPoints = 1;
  for (int i = 1; i < l; i++)
     {int j = 0;
      while (j < fs[fr].nPoints)
	switch (dominates2way(fs[fr].points[i], fs[fr].points[j], n-2))
	  {case  0: j++; break;
	   case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
	            // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
		   {
	            t = fs[fr].points[j];
		    fs[fr].points[j] = fs[fr].points[i]; 
		    fs[fr].points[i] = t; 
		    while(j < fs[fr].nPoints - 1 && dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].nPoints - 1], n-1))
		      fs[fr].nPoints--;
		    int k = j+1; 
		    while (k < fs[fr].nPoints)
		      if(dominates1way(fs[fr].points[j], fs[fr].points[k], n-2))
			{t = fs[fr].points[k];
			 fs[fr].nPoints--;
			 fs[fr].points[k] = fs[fr].points[fs[fr].nPoints]; 
			 fs[fr].points[fs[fr].nPoints] = t; 
			}
		      else
			k++;
		   }
	   default: j = fs[fr].nPoints + 1;
	  }
      if (j == fs[fr].nPoints) 
	{t = fs[fr].points[fs[fr].nPoints]; 
	 fs[fr].points[fs[fr].nPoints] = fs[fr].points[i]; 
	 fs[fr].points[i] = t; 
	 fs[fr].nPoints++;}
     }
  safe = WORSE(l,fs[fr].nPoints);
  for (int i = l; i < p; i++)
     {int j = 0;
      while (j < safe)
	if(dominates1way(fs[fr].points[j], fs[fr].points[i], n-2))
	  j = fs[fr].nPoints + 1;
	else
	  j++;
      while (j < fs[fr].nPoints)
	switch (dominates2way(fs[fr].points[i], fs[fr].points[j], n-1))
	  {case  0: j++; break;
	   case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
	            // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
		   {
	            t = fs[fr].points[j];
		    fs[fr].points[j] = fs[fr].points[i]; 
		    fs[fr].points[i] = t; 
		    while(j < fs[fr].nPoints - 1 && dominates1way(fs[fr].points[j], fs[fr].points[fs[fr].nPoints - 1], n-1))
		      fs[fr].nPoints--;
		    int k = j+1; 
		    while (k < fs[fr].nPoints)
		      if(dominates1way(fs[fr].points[j], fs[fr].points[k], n-1))
			{t = fs[fr].points[k];
			 fs[fr].nPoints--;
			 fs[fr].points[k] = fs[fr].points[fs[fr].nPoints]; 
			 fs[fr].points[fs[fr].nPoints] = t; 
			}
		      else
			k++;
		   }
	   default: j = fs[fr].nPoints + 1;
	  }
      if (j == fs[fr].nPoints) 
	{t = fs[fr].points[fs[fr].nPoints]; 
	  fs[fr].points[fs[fr].nPoints] = fs[fr].points[i]; 
	  fs[fr].points[i] = t; 
	  fs[fr].nPoints++;}
     }
  fr++;
}


double hv2(FRONT ps, int k)
// returns the hypervolume of ps[0 .. k-1] in 2D 
// assumes that ps is sorted improving
{
  double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1]; 
  for (int i = 1; i < k; i++) 
    volume += ps.points[i].objectives[1] * 
              (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
  return volume;
}


double inclhv(POINT p)
// returns the inclusive hypervolume of p
{
  double volume = 1;
  for (int i = 0; i < n; i++) 
    volume *= p.objectives[i];
  return volume;
}


double inclhv2(POINT p, POINT q)
// returns the hypervolume of {p, q}
{
  double vp  = 1; double vq  = 1;
  double vpq = 1;
  for (int i = 0; i < n; i++) 
    {
      vp  *= p.objectives[i];
      vq  *= q.objectives[i];
      vpq *= WORSE(p.objectives[i],q.objectives[i]);
    }
  return vp + vq - vpq;
}


double inclhv3(POINT p, POINT q, POINT r)
// returns the hypervolume of {p, q, r}
{
  double vp   = 1; double vq   = 1; double vr   = 1;
  double vpq  = 1; double vpr  = 1; double vqr  = 1;
  double vpqr = 1;
  for (int i = 0; i < n; i++) 
    {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
	if (BEATS(q.objectives[i],r.objectives[i]))
	{
	  vpq  *= q.objectives[i];
	  vpr  *= r.objectives[i];
	  vqr  *= r.objectives[i];
	  vpqr *= r.objectives[i];
	}
	else
	{
	  vpq  *= q.objectives[i];
	  vpr  *= WORSE(p.objectives[i],r.objectives[i]);
	  vqr  *= q.objectives[i];
	  vpqr *= q.objectives[i];
	}
      else
	if (BEATS(p.objectives[i],r.objectives[i]))
	{
	  vpq  *= p.objectives[i];
	  vpr  *= r.objectives[i];
	  vqr  *= r.objectives[i];
	  vpqr *= r.objectives[i];
	}
	else
	{
	  vpq  *= p.objectives[i];
	  vpr  *= p.objectives[i];
	  vqr  *= WORSE(q.objectives[i],r.objectives[i]);
	  vpqr *= p.objectives[i];
	}
    }
  return vp + vq + vr - vpq - vpr - vqr + vpqr;
}


double inclhv4(POINT p, POINT q, POINT r, POINT s)
// returns the hypervolume of {p, q, r, s}
{
  double vp    = 1; double vq   = 1; double vr   = 1; double vs   = 1;
  double vpq   = 1; double vpr  = 1; double vps  = 1; double vqr  = 1; double vqs = 1; double vrs = 1; 
  double vpqr  = 1; double vpqs = 1; double vprs = 1; double vqrs = 1; 
  double vpqrs = 1; 
  for (int i = 0; i < n; i++) 
    {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      vs *= s.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
	if (BEATS(q.objectives[i],r.objectives[i]))
	  if (BEATS(r.objectives[i],s.objectives[i]))
	    {
	      vpq *= q.objectives[i];
	      vpr *= r.objectives[i];
	      vps *= s.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= r.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
	      vpq *= q.objectives[i];
	      vpr *= r.objectives[i];
	      vps *= WORSE(p.objectives[i],s.objectives[i]);
	      vqr *= r.objectives[i];
	      vqs *= z1;
	      vrs *= r.objectives[i];
	      vpqr *= r.objectives[i];
	      vpqs *= z1;
	      vprs *= r.objectives[i];
	      vqrs *= r.objectives[i];
	      vpqrs *= r.objectives[i];
	    }
	else
	  if (BEATS(q.objectives[i],s.objectives[i]))
	    {
	      vpq *= q.objectives[i];
	      vpr *= WORSE(p.objectives[i],r.objectives[i]);
	      vps *= s.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= q.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      vpq *= q.objectives[i];
	      vpr *= z1;
	      vps *= WORSE(p.objectives[i],s.objectives[i]);
	      vqr *= q.objectives[i];
	      vqs *= q.objectives[i];
	      vrs *= WORSE(r.objectives[i],s.objectives[i]);
	      vpqr *= q.objectives[i];
	      vpqs *= q.objectives[i];
	      vprs *= WORSE(z1,s.objectives[i]);
	      vqrs *= q.objectives[i];
	      vpqrs *= q.objectives[i];
	    }
      else
	if (BEATS(q.objectives[i],r.objectives[i]))
	  if (BEATS(p.objectives[i],s.objectives[i]))
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= z1;
	      vps *= s.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= z2;
	      vpqr *= z1;
	      vpqs *= s.objectives[i];
	      vprs *= z2;
	      vqrs *= z2;
	      vpqrs *= z2;
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= z1;
	      vps *= p.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= WORSE(q.objectives[i],s.objectives[i]);
	      vrs *= z2;
	      vpqr *= z1;
	      vpqs *= p.objectives[i];
	      vprs *= z1;
	      vqrs *= z2;
	      vpqrs *= z1;
	    }
	else
	  if (BEATS(p.objectives[i],s.objectives[i]))
	    {
	      vpq *= p.objectives[i];
	      vpr *= p.objectives[i];
	      vps *= s.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= p.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= p.objectives[i];
	      vps *= p.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= z1;
	      vrs *= WORSE(r.objectives[i],s.objectives[i]);
	      vpqr *= p.objectives[i];
	      vpqs *= p.objectives[i];
	      vprs *= p.objectives[i];
	      vqrs *= z1;
	      vpqrs *= p.objectives[i];
	    }
    }
  return vp + vq + vr + vs - vpq - vpr - vps - vqr - vqs - vrs + vpqr + vpqs + vprs + vqrs - vpqrs;
}


double exclhv(FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[0 .. p-1] 
{
  makeDominatedBit(ps, p);
  double volume = inclhv(ps.points[p]) - hv(fs[fr - 1]);
  fr--;
  return volume;
}


double hv(FRONT ps)
// returns the hypervolume of ps[0 ..] 
{
  // process small fronts with the IEA 
  switch (ps.nPoints)
    {case 1: return inclhv (ps.points[0]); 
     case 2: return inclhv2(ps.points[0], ps.points[1]); 
     case 3: return inclhv3(ps.points[0], ps.points[1], ps.points[2]); 
     case 4: return inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]); 
    }

  // these points need sorting 
  qsort(&ps.points[safe], ps.nPoints - safe, sizeof(POINT), greater); 
  // n = 2 implies that safe = 0 
  if (n == 2) return hv2(ps, ps.nPoints); 
  // these points don't NEED sorting, but it helps 
  qsort(ps.points, safe, sizeof(POINT), greaterabbrev); 

  if (n == 3 && safe > 0) 
    {
      double volume = ps.points[0].objectives[2] * hv2(ps, safe); 
      n--;
      for (int i = safe; i < ps.nPoints; i++)
	// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
		volume += ps.points[i].objectives[n] * exclhv(ps, i);
      n++; 
      return volume;
    }
  else
    {
      double volume = inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]); 
      n--;
      for (int i = 4; i < ps.nPoints; i++)
	// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
		volume += ps.points[i].objectives[n] * exclhv(ps, i);
      n++; 
      return volume;
    }
}


int main() 
{
  dim = 10;	//number of objectives
  printf("obj=%d\n",dim);
  	// initialise the reference point
  ref.objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * dim);
  //double rep[5] = {6092, 6019, 6629, 5763, 4090};
  for(int i=0;i<dim;++i)
  {
  	ref.objectives[i]=2.2*(i+1);//rep[i];
  }    
  double refsum=1.0;
  for (int i = 0; i < dim; i++) 
  {
  	refsum*=ref.objectives[i];
  }
  
  //FILECONTENTS *f = readFile("E:\\matlab\\MODEA-test results\\2D-3D\\UF10");
  //FILECONTENTS *f = readFile("E:\\matlab\\MOEA-D-test results\\MOEA-D-DE\\2D-3D\\UF10");
  //FILECONTENTS *f = readFile("E:\\matlab\\CMPSO-test results\\2D-3D\\UF10");
  //FILECONTENTS *f = readFile("E:\\matlab\\CMODE\\CMODE\\2D-3D\\UF10");
  
  //FILECONTENTS *f = readFile("E:\\matlab\\NSGA3-test results\\5D-8D-10D(140w)\\UF13");
  //FILECONTENTS *f = readFile("E:\\matlab\\MOEA-D-test results\\MOEA-D-DE\\5D-8D-10D(140w)\\UF13");
 // FILECONTENTS *f = readFile("E:\\matlab\\SPEA2+SDE-test results\\5D-8D-10D(140w)\\UF13");
 //FILECONTENTS *f=readFile("./obj10/WFG3/chebyshev_angle_share_finalfit_g=300.txt");
 FILECONTENTS *f=readFile("obj10/WFG1/origin_finalfit_g=300.txt");
 //FILECONTENTS *f = readFile("E:\\matlab\\CMODE\\CMODE\\5D-8D-10D(140w)\\UF13");
  
  //FILECONTENTS *f = readFile("E:\\matlab\\CMDE+KNN-SDE-test results\\5D-8D-10D\\UF13");
  // find the biggest fronts
  for (int i = 0; i < f->nFronts; i++)
    {if (f->fronts[i].nPoints > maxm) {maxm = f->fronts[i].nPoints;}
     if (f->fronts[i].n       > maxn) {maxn = f->fronts[i].n;}
     //printf("%d %d\n", maxm, maxn);
    }
  // allocate memory
  int maxdepth = maxn - 2; 
  fs = (FRONT*)malloc(sizeof(FRONT) * maxdepth);
  for (int i = 0; i < maxdepth; i++) 
    {
    	fs[i].points = (POINT*)malloc(sizeof(POINT) * maxm); 
     	for (int j = 0; j < maxm; j++) 
     	{
       		fs[i].points[j].objectives = (OBJECTIVE*)malloc(sizeof(OBJECTIVE) * (maxn - i - 1));
   		}
    }
  /*
  if (argc == 2)
    {printf("No reference point provided: using the origin\n");
     for (int i = 0; i < maxn; i++) ref.objectives[i] = 0;
    }
  else if (argc - 2 != maxn)
    {printf("Your reference point should have %d values\n", maxn);
     return 0;
    }
  else 
  for (int i = 2; i < argc; i++) ref.objectives[i - 2] = atof(argv[i]);
  */

  // modify the objective values relative to the reference point 
  double startTime=clock();
  for (int i = 0; i < f->nFronts; i++) 
  {
    for(int j = 0; j < f->fronts[i].nPoints; j++)
    {
      for(int k = 0; k < f->fronts[i].n; k++)
      {
	      f->fronts[i].points[j].objectives[k] = fabs(f->fronts[i].points[j].objectives[k] - ref.objectives[k]);
	  }
	}
  }
  
  int tim=f->nFronts;
  double frontHv[tim];
  double maxhv=-1;
  double minhv=2;
  double averhv=0;
  for (int i = 0; i < tim; i++) 
    {      
 
      n = f->fronts[i].n;
      safe = 0;
      frontHv[i]=hv(f->fronts[i])/refsum;
      printf("hv(%d) = %.10f\n", i+1, frontHv[i]); 
      
      averhv+=frontHv[i];
      if(frontHv[i]>maxhv){maxhv=frontHv[i];}
      if(frontHv[i]<minhv){minhv=frontHv[i];}
    }
    averhv/=tim;
    double stdhv=0;
    for(int i=0;i<tim;++i)
    {
    	stdhv+=pow(frontHv[i]-averhv, 2.0);
    }
    stdhv=sqrt(stdhv/tim);
    
    double endTime=clock();
    double runTime=(endTime-startTime)/CLOCKS_PER_SEC;
    printf("\nrun time is %f(s)\n\n", runTime);
    printf("maxHV:  %.6lf\n",maxhv);
    printf("averHV: %.6lf\n",averhv);
    printf("minHV:  %.6lf\n",minhv);
    printf("stdHV:  %.6lf\n",stdhv);
    
    system("pause");
  return 0;
}

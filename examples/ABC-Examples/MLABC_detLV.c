/* LIBABC: approximate Bayesian Computation
 * Copyright (C) 2016  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mcl.h"
#include <time.h>

/* Example of ABC-PCR samples for a deteministic Lotka-Volterra system from Toni et al.*/
ABC_Parameters abc_p;
MLMC_Parameters mlmc_p;


Dataset data;

/*generate grid for indicator thresholds*/
ndgrid GenNdGrid(int,double *,double *,int,double *);

/*discrepency metric */
double rho(Dataset *data, Dataset *data_s)
{
    double *Y,*Y_s;
    double d;
    unsigned int nt,n,i,j;
    Y_s = (double *)data_s->fields[0].data_array;
    Y = (double *)data->fields[0].data_array;

    nt = (unsigned int)data->fields[0].numRows;
    n = (unsigned int)data->fields[0].numCols;

    d = 0;

    for (j=0;j<n;j++)
    {
        for (i=0;i<nt;i++)
        {
            d += (Y[j*nt + i] - Y_s[j*nt + i])*(Y[j*nt + i] - Y_s[j*nt + i]);
        }
    }
    return d;
}

/*prior sampler a ~ U(-10,10) b ~ U(-10,10)*/
int prior(unsigned int dim, unsigned int nsamples, double* support,double * theta)
{
    unsigned int i,j;
    if (support == NULL)
    {
        for (i=0;i<nsamples;i++)
        {
            for (j=0;j<dim;j++)
            {
                theta[i*dim + j] = durngus(0.5,1.2); 
            }
        }
    }
    else
    {
        for (i=0;i<nsamples;i++)
        {
            for (j=0;j<dim;j++)
            {
                theta[i*dim + j] = durngus(support[j],support[j+dim]); 
            }
        }
 
    }
    return 0;
}

/* Deterministic Lotka-Volterra preditor/prey model
 *
 */
void detLV(double *Y, unsigned int n, double *params,unsigned int m, double t,double *f_r)
{
   double a,b,x,y;
   a = params[0];
   b = params[1];
   x = Y[0];
   y = Y[1];

   f_r[0] = a*x - x*y;
   f_r[1] = b*x*y - y;
}

#define TSTEP 1.5e-3
/*simulation*/
int simulate(void *sim,double * theta, Dataset * data_s)
{
    unsigned int n,m,nt,i;
    double *Y_r;
    double Y0[2] = {1.0,0.5};
    double T[8] = {1.1,2.4,3.9,5.6,7.5,9.6,11.9,14.4}; 
    unsigned int d[2] = {0,1};
    double h;
    h = TSTEP;

    /*intial conditions*/
    Y_r = (double*)data_s->fields[0].data_array;
    nt = data_s->fields[0].numRows;
    n = data_s->fields[0].numCols;

    if (n != 2 || nt != 8)
    {
        fprintf(stderr,"Fatal Error: Invalid dataset [n,nt = %d %d].\n",n,nt);
        exit(1);
    }

    m = 2; /*params are a, b */
    
    /*solve ODE with RK4 method*/
    drk4s(m,n,nt,T,theta,Y0,&detLV,n,d,h,Y_r);

    sim_counter++;
    return 0;
}

double deltas[2];

/*smoothing function*/
SSAL_real_t g(int d, SSAL_real_t *y, SSAL_real_t *x)
{
    int j;
    SSAL_real_t res = 1;
    SSAL_real_t s[255];

    for (j=0;j<d;j++)
    {
        s[j] = (y[j] -x[j])/deltas[j];    
        //s[j] = (y[j] -x[j]);    
        /*based on Giles method*/  
        res *= (s[j] <= 1) ? ((s[j] >= -1) ? 0.625*s[j]*s[j]*s[j] - 1.125*s[j] + 0.5 : 1 ) : 0;
    }
    return res;
}


int main(int argc, char ** argv)
{
    double *Nl_scale;
    double *theta,*weights;
    CDF_estimate M;
    clock_t start_t, end_t;
    double time;
    unsigned int i,l;
    double *Y_d;
    /*only 2 args*/
    if (argc  < 2)
    {
        fprintf(stderr,"Usage: %s N\n",argv[0]);
    }
    else
    {
        /*set up ABC params*/
        abc_p.nacc = (unsigned int)atoi(argv[1]);
        abc_p.eps = 4.3;
        abc_p.nmax = 0;
        abc_p.k = 2;
        abc_p.support = (double*)malloc(4*sizeof(double));
        abc_p.support[0] = 0.5;
        abc_p.support[1] = 0.5;
        abc_p.support[2] = 1.2;
        abc_p.support[3] = 1.2;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
 
    /*set up MLMC parameters*/
    mlmc_p.L = 4;
    mlmc_p.eps_l = (double*)malloc(5*sizeof(double));
    mlmc_p.eps_l[0] = 30.0;
    mlmc_p.eps_l[1] = 16.0; 
    mlmc_p.eps_l[2] = 6.0;
    mlmc_p.eps_l[3] = 5.0;
    mlmc_p.eps_l[4] = 4.3;
    mlmc_p.Nl = (unsigned int*)malloc((mlmc_p.L+1)*sizeof(unsigned int));
    Nl_scale = (double*)malloc((mlmc_p.L+1)*sizeof(double));
    mlmc_p.target_RMSE = 0.1;
    mlmc_p.presample = 1; 
    mlmc_p.presample_trials = 100; 
    
    /*initialise grid: we define it only in */
    //double lb[2] = {0.8,0.8};
    //double ub[2] = {1.15,1.15};
    //deltas[0] = (ub[0] - lb[0])/999.0; 
    //deltas[1] = (ub[1] - lb[1])/999.0; 
    //M.G = GenNdGrid(2,lb,ub,1000,deltas);
    deltas[0] = (abc_p.support[2] - abc_p.support[0])/999.0; 
    deltas[1] = (abc_p.support[3] - abc_p.support[1])/999.0; 
   // deltas[1] = (ub[1] - lb[1])/999.0; 
    M.G = GenNdGrid(abc_p.k,abc_p.support,abc_p.support + abc_p.k,1000,deltas);
    M.F = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.V = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.g = &g; 


    
    /*initialise our RNG library*/
    SSAL_Initialise(argc,argv);

    /*build data template*/
    data.numFields = 1 ;
    data.fields = (field *)malloc(sizeof(field));
    data.fields[0].numCols = 2;
    data.fields[0].numRows = 8;
    data.fields[0].type = REAL64_DATA;
    data.fields[0].data_array = malloc(16*sizeof(double));
    data.fields[0].numBytes = 16*sizeof(double);
   
    /* populate data (pre-computed)*/
    Y_d = (double*)data.fields[0].data_array;
    Y_d[0] = 1.87; Y_d[8] = 0.49;
    Y_d[1] = 0.65; Y_d[9] = 2.62;
    Y_d[2] = 0.22; Y_d[10] = 1.54;
    Y_d[3] = 0.31; Y_d[11] = 0.02;
    Y_d[4] = 1.64; Y_d[12] = 1.14;
    Y_d[5] = 1.15; Y_d[13] = 1.68;
    Y_d[6] = 0.24; Y_d[14] = 1.07;
    Y_d[7] = 2.91; Y_d[15] = 0.88;

    /*allocate output array*/
    sim_counter = 0; /*for performance metric*/
    /*determine sample number scaling for MLMC*/
    start_t = clock();
    dmlabcnls(abc_p,mlmc_p.presample_trials,mlmc_p,&data,&M,mlmc_p.Nl,Nl_scale);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stderr,"\"Level\",\"N_l\",\"N_l_s\",\"SEC\",\"NSIMS\"\n");
    for (l=0;l<=mlmc_p.L;l++)
    {
        fprintf(stderr,"%u,%u,%g,%g,%u\n",l,mlmc_p.Nl[l],Nl_scale[l],time,sim_counter);
    }
    for (l=0;l<=mlmc_p.L;l++)
    {
        mlmc_p.Nl[l] = (Nl_scale[l]/Nl_scale[mlmc_p.L])*abc_p.nacc;
    }
    fprintf(stderr,"\"Level\",\"N_l\",\"N_l_s\",\"SEC\",\"NSIMS\"\n");
    for (l=0;l<=mlmc_p.L;l++)
    {
        fprintf(stderr,"%u,%u,%g,%g,%u\n",l,mlmc_p.Nl[l],Nl_scale[l],time,sim_counter);
    }
    /*run MLABC*/
    sim_counter = 0; /*for performance metric*/
    start_t = clock();
    dmlabccdfs(abc_p,mlmc_p,&data,&M);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write ouput*/
    fprintf(stdout,"\"a\",\"b\",\"F\",\"SEC\",\"NSIMS\"\n");
    for (i=0;i<M.G.numPoints;i++)
    {
        fprintf(stdout,"%g,%g,%g,%g,%u\n",M.G.coords[i*2],M.G.coords[i*2+1],M.F[i],time,sim_counter);
    }

    
    return 0;
}

/*generate a regular N-D grid*/
ndgrid GenNdGrid(int d, SSAL_real_t *S_l,SSAL_real_t *S_u,int D,double * deltas)
{
    SSAL_real_t * axes;
    size_t *index;
    int i,j;
    ndgrid G; 

    index = (size_t *)malloc(d*sizeof(size_t));
    axes = (SSAL_real_t *)malloc(d*D*sizeof(SSAL_real_t));
    G.dim = d;
    G.D = D;
   
    memcpy(G.deltas,deltas,G.dim*sizeof(double));

    for (j=0;j<d;j++)
    {
        for (i=0;i<D;i++)
        {
            axes[j*D+i] = (((SSAL_real_t)i)*(S_u[j] - S_l[j]))/((SSAL_real_t)(D-1)) + S_l[j];
        }
    }

    G.numPoints = 1;
    for (j=0;j<d;j++)
    {
        G.numPoints *= D;
    }

    G.coords = (SSAL_real_t *)malloc(G.numPoints*d*sizeof(SSAL_real_t));
    /*init index*/
    for (j=0;j<d;j++)
    {
        index[j] = 0;
        G.coords[j] = axes[j*D + index[j]];
    }
    
    int c; /*carry flag*/
    /*generate grid points*/
    for (i=1;i<G.numPoints;i++)
    {
        c = 1;
        for (j=0;j<d;j++) /*base D adder with carry*/
        {
            index[j] = (index[j] + c) % D;
            c = (index[j] == 0 && c != 0);
        }
        /* generate point*/

        for (j=0;j<d;j++)
        {
            G.coords[i*d + j] = axes[j*D + index[j]];
        }
    }

    /*compute offsets and increments (for marginals) */
    for (j=0;j<d;j++)
    {
        int mult;
        mult = 1;
        G.offsets[j] = 0;
        for (i=0;i<d;i++)
        {
            if (i == j)
            {
                G.incs[j] = mult;
            }
            else
            {
                G.offsets[j] += (D-1)*mult;
            }
            mult *= D;
        }
    }

    free(index);
    free(axes);
    return G;
    
}



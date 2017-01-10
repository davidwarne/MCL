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
SMC_Parameters smc_p;


Dataset data;

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
    for (i=0;i<nsamples;i++)
    {
        for (j=0;j<dim;j++)
        {
           theta[i*dim + j] = durngus(0.5,1.2); 
        }
    }
    return 0;
}

/* for SMC we also need densities*/
double priorPDF(unsigned int dim,double *theta)
{
    unsigned int i;
    double p;
    p = 1;
    for (i=0;i<dim;i++)
    {
        p *= (1.0/0.7);
    }
    return p; 
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

#define SIGMA 0.1
/*transition kernel and PDF*/
int kernel(unsigned int dim,double *theta, double *theta_prop)
{
    unsigned int i;
    for (i=0;i<dim;i++)
    {
        theta_prop[i] = theta[i] + SIGMA*durngus(-1.0,1.0);
    }
    return 0;
}
double kernelPDF(unsigned int dim,double *theta, double *theta_prop)
{
    unsigned int i;
    double p;
    p = 1;
    for (i=0;i<dim;i++)
    {
        p *= (1.0/(2.0*SIGMA));
    }
    return p; 
}

int main(int argc, char ** argv)
{
    double *theta,*weights;
    clock_t start_t, end_t;
    double time;
    unsigned int i;
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
        abc_p.support = NULL;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = &priorPDF;
        abc_p.s = &simulate;
    }
 
    /*set up SMC parameters*/
    smc_p.T = 5;
    smc_p.eps_t = (double*)malloc(5*sizeof(double));
    smc_p.eps_t[0] = 30.0;
    smc_p.eps_t[1] = 16.0; 
    smc_p.eps_t[2] = 6.0;
    smc_p.eps_t[3] = 5.0;
    smc_p.eps_t[4] = 4.3;
    smc_p.E = ((double)abc_p.nacc)/2.0;
    smc_p.q = &kernel;
    smc_p.qd = &kernelPDF;
    
    
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
    theta = (double *)malloc(abc_p.nacc*abc_p.k*sizeof(double));
    weights = (double *)malloc(abc_p.nacc*sizeof(double));

    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC-MCMC*/
    start_t = clock();
    dabcpcr(abc_p,smc_p,&data,theta,weights,NULL);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stdout,"\"Particle\",\"a\",\"b\",\"weights\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%g,%u,%u\n",i,theta[i*2],theta[i*2+1],weights[i],time,abc_p.nacc,sim_counter);
    }

    return 0;
}


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

#include "libabc.h"
#include <time.h>

/* Example of ABC-PCR samples for toy problem from Sisson et al.*/
ABC_Parameters abc_p;
SMC_Parameters smc_p;


Dataset data; /*not really used except as a template for simulated data*/
/*discrepency metric (data implicitly used)*/
double rho(Dataset *data, Dataset *data_s)
{
    double *x;
    x = (double *)data_s->fields[0].data_array;

    if (durngus(0.0,1.0) < 0.5)
    {
        return fabs(x[0]); 
    }
    else
    {
        size_t n;
        double m;
        int i;
        n = data_s->fields[0].numCols; 
        m = 0.0;
        for (i=0;i<n;i++)
        {
            m += x[i];
        }
        m /= (double)n;
        return fabs(m);
    }
}

/*prior sampler*/
int prior(unsigned int dim, unsigned int nsamples, double* support,double * theta)
{
    unsigned int i;
    for (i=0;i<nsamples;i++)
    {
        theta[i] = durngus(-10.0,10.0);
    }
    return 0;
}

/* for MCMC we also need densities*/
double priorPDF(unsigned int dim,double *theta)
{
    return 1.0/20.0;
}

/*simulation*/
int simulate(void *sim,double * theta, Dataset * data_s)
{
    unsigned int n,i;
    double *x;
    x = (double*)data_s->fields[0].data_array;
    n = data_s->fields[0].numCols;
    
    for (i=0;i<n;i++)
    {
        x[i] = durngns(theta[0],1.0);
        //printf("X[%d] ~ N(%f,1) = %f\n",i,theta[0],x[i]);
    }
    sim_counter++;
    return 0;
}

#define SIGMA 0.15
/*transition kernel and PDF*/
int kernel(unsigned int dim,double *theta, double *theta_prop)
{
    theta_prop[0] = durngns(theta[0],SIGMA);
}
double kernelPDF(unsigned int dim,double *theta, double *theta_prop)
{
    double a,b;
    a = 1.0/(SIGMA*sqrt(2.0*M_PI));
    b = ((theta_prop[0] - theta[0])/SIGMA);
    return a*exp(-0.5*b*b);
}

int main(int argc, char ** argv)
{
    double *theta,*weights;
    clock_t start_t, end_t;
    double time;
    unsigned int i;
    double theta0;
    theta0 =0;
    /*only 2 args*/
    if (argc  < 3)
    {
        fprintf(stderr,"Usage: %s N eps\n",argv[0]);
    }
    else
    {
        /*set up ABC params*/
        abc_p.nacc = (unsigned int)atoi(argv[1]);
        abc_p.eps = (double)atof(argv[2]);
        abc_p.nmax = 0;
        abc_p.k = 1;;
        abc_p.support = NULL;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = &priorPDF;
        abc_p.s = &simulate;
    }
 
    /*set up SMC parameters*/
    smc_p.T = 3 ;
    smc_p.eps_t = (double*)malloc(3*sizeof(double));
    smc_p.eps_t[0] = 2.0;
    smc_p.eps_t[1] = 0.5; 
    smc_p.eps_t[2] = 0.025;
    smc_p.E = ((double)abc_p.nacc)/2.0;
    smc_p.q = &kernel;
    smc_p.qd = &kernelPDF;
    
    
    /*initialise our RNG library*/
    SSAL_Initialise(argc,argv);

    /*build data template*/
    data.numFields = 1 ;
    data.fields = (field *)malloc(sizeof(field));
    data.fields[0].numCols = 100;
    data.fields[0].numRows = 1;
    data.fields[0].type = REAL64_DATA;
    data.fields[0].data_array = malloc(100*sizeof(double));
    data.fields[0].numBytes = 100*sizeof(double);
    
    /*allocate output array*/
    theta = (double *)malloc(abc_p.nacc*sizeof(double));
    weights = (double *)malloc(abc_p.nacc*sizeof(double));

    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC-MCMC*/
    start_t = clock();
    dabcpcr(abc_p,smc_p,&data,theta,weights,NULL);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stdout,"\"Particle\",\"theta\",\"weights\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%u,%u\n",i,theta[i],weights[i],time,abc_p.nacc,sim_counter);
    }

    return 0;
}


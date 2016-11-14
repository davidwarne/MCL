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

/* Example of ABC samples for toy problem from Sisson et al.*/
ABC_Parameters abc_p;

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
    }
    sim_counter++;
    return 0;
}

int main(int argc, char ** argv)
{
    double *theta;
    clock_t start_t, end_t;
    double time;
    unsigned int i;

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
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
    
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

    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC rejection*/
    start_t = clock();
    dabcrs(abc_p,&data,theta,NULL);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stdout,"\"Realisation\",\"theta\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%u,%u\n",i,theta[i],time,abc_p.nacc,sim_counter);
    }

    return 0;
}


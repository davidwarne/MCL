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

/* Example of ABC-rejection samples for TB transmission model of Tanaka et al.*/
ABC_Parameters abc_p;

Dataset data; /*not really used except as a template for simulated data*/
/*discrepency metric (data implicitly used)*/
double rho(Dataset *data, Dataset *data_s)
{
    double *gH_s,*gH_d;
    double n,r;
    gH_s = (double *)data_s->fields[0].data_array;
    gH_d = (double *)data->fields[0].data_array;

    n = gH_d[2];

    r = (1.0/n)*fabs(gH_s[0] - gH_d[0]) + fabs(gH_s[1]-gH_s[1]);
    return r;
}
/*prior sampler*/
int prior(unsigned int dim, unsigned int nsamples, double* support,double * theta)
{
    theta[0] = durngus(0.0,5.0);/*alpha: birth rate*/
    theta[1] = durngus(0.0,theta[0]);/*delta: death rate*/
    theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
    while (theta[2] < 0.0)
    {
        theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
    }
    return 0;
}

#define BIRTH 0
#define DEATH 1
#define MUTATE 2
/*simulation*/
int simulate(void *sim,double * theta, Dataset * data_s)
{
    /*probably set these to get from sim parameters*/
    unsigned int N_stop;
    unsigned int n;
    double *X; /*X[i] number of infections from TB bacterium of genotype i */
    double *x; /* sub-population sample of X */
    unsigned int G; /*number of genotypes of the TB bacterium*/
    double N; /* total infection cases*/
    double sum_theta;
    unsigned int i,j,event;
    double *gH;
   
    n = 473;
    N_stop = 10000;

    gH = (double*)data_s->fields[0].data_array;

    /*initialise  model*/
    X = (double *)malloc(N_stop*sizeof(double));
    x = (double *)malloc(N_stop*sizeof(double));
    memset(X,0,N_stop*sizeof(double));
    memset(x,0,N_stop*sizeof(double));
    /*start with a single case */
    X[0] = 1;
    N = 1;
    G = 1;

    sum_theta = theta[0] + theta[1] + theta[2];/*alpha + delta + theta*/
    
    /*Gillespie simulations*/
    while (N > 0.0 && N < (double)N_stop)
    {
        /*select event */
        event = durngpmfs(3,theta,sum_theta);

        /*select genotype for event*/
        i = durngpmfs(G,X,N);

        /*simulate the next event*/
        switch (event)
        {
            case BIRTH:
                X[i] += 1.0;
                N += 1.0;
                break;
            case DEATH:
                X[i] -= 1.0;
                N -= 1.0;
                if (X[i] == 0.0)
                {
                    for (j=i;j<G-1;j++)
                    {
                        X[j] = X[j+1];
                    }
                    G--;
                }
                break;
            case MUTATE:
                X[i] -= 1.0;
                if (X[i] == 0.0)
                {
                    for (j=i;j<G-1;j++)
                    {
                        X[j] = X[j+1];
                    }
                    G--;
                }

                X[G] = 1.0;
                G++;
                break;
        }
    }

    /*population died out */
    if (N <= 0.0)
    {
        gH[0] = 0;
        gH[1] = 0;
        free(X);
        free(x);
        return 0;
    }

    for (i=0;i<n;i++)
    {
        j = durngpmfs(G,X,N);
        x[j]+= 1.0;
        X[j]-=1.0;
        N -= 1.0;
    }

    /*compute the number of distinct genotypes and genetic diversity*/
    gH[0] = 0.0;
    gH[1] = 1.0;
    for (i=0;i<G;i++)
    {
        if (x[i] > 0.0)
        {
            gH[0] += 1.0;
            gH[1] -= (x[i]*x[i])/((double) n*n);
        }
    }
    gH[2] = (double)n;
    sim_counter++;
    free(X);
    free(x);
    return 0;
}

int main(int argc, char ** argv)
{
    double *theta;
    double * gH_d;
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
        abc_p.k = 3;;
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
    data.fields[0].numCols = 3;
    data.fields[0].numRows = 1;
    data.fields[0].type = REAL64_DATA;
    data.fields[0].data_array = malloc(3*sizeof(double));
    data.fields[0].numBytes = 3*sizeof(double);
    /*assign data*/
    gH_d = (double *)data.fields[0].data_array;
    /* IS6110 fingerprint from Small et al.*/
    /*30^1 23^1 15^1 10^1 8^1 5^2 4^4 3^13 2^20 1^282*/
    gH_d[0] = 326.0;
    gH_d[1] = 0.9892236;
    gH_d[2] = 473.0;


    /*allocate output array*/
    theta = (double *)malloc(abc_p.nacc*3*sizeof(double));

    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC rejection*/
    start_t = clock();
    dabcrs(abc_p,&data,theta,NULL);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    
    /*write output*/
    fprintf(stdout,"\"Interations\",\"alpha\",\"delta\",\"theta\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%g,%u,%u\n",i,theta[i*3],theta[i*3 + 1],theta[i*3 + 2],time,abc_p.nacc,sim_counter);
    }

    free(theta);
    return 0;
}


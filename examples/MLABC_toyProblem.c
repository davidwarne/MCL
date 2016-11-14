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
MLMC_Parameters mlmc_p;


Dataset data; /*not really used except as a template for simulated data*/
/*generate grid for indicator thresholds*/
ndgrid GenNdGrid(int,double *,double *,int,double *);
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
    if (support == NULL)
    {
        for (i=0;i<nsamples;i++)
        {
            theta[i] = durngus(-10.0,10.0);
        }
    }
    else
    {
        for (i=0;i<nsamples;i++)
        {
            theta[i] = durngus(support[0],support[1]);
        }
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
        //printf("X[%d] ~ N(%f,1) = %f\n",i,theta[0],x[i]);
    }
    sim_counter++;
    return 0;
}

double deltas[1];

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
    CDF_estimate M;
    clock_t start_t, end_t;
    double time;
    unsigned int i,l;
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
        abc_p.k = 1;
        abc_p.support = (double*)malloc(2*sizeof(double));
        abc_p.support[0] = -10.0;
        abc_p.support[1] = 10.0;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
 
    /*set up MLMC parameters*/
    mlmc_p.L = 2 ;
    mlmc_p.eps0 = 2.0;
    mlmc_p.eps_l = (double*)malloc(3*sizeof(double));
    mlmc_p.Nl = (unsigned int*)malloc(3*sizeof(unsigned int));
    Nl_scale = (double*)malloc(3*sizeof(double));
    mlmc_p.eps_l[0] = 2.0;
    mlmc_p.eps_l[1] = 0.5; 
    mlmc_p.eps_l[2] = 0.025;
    mlmc_p.target_RMSE = 0.1;
    mlmc_p.presample = 1; 
    mlmc_p.presample_trials = 100; 

    /*initialise grid*/
    deltas[0] = (abc_p.support[1] - abc_p.support[0])/999.0; 
    M.G = GenNdGrid(1,abc_p.support,abc_p.support+1,1000,deltas);
    M.F = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.V = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.g = &g; 

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

    sim_counter = 0; /*for performance metric*/
    /*determine sample number scaling for MLMC*/
    start_t = clock();
    dmlabcnls(abc_p,mlmc_p.presample_trials,mlmc_p,&data,&M,mlmc_p.Nl,Nl_scale);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    ///*write output*/
    fprintf(stderr,"\"Level\",\"N_l\",\"N_l_s\",\"SEC\",\"NSIMS\"\n");
    for (l=0;l<=mlmc_p.L;l++)
    {
        fprintf(stderr,"%u,%u,%g,%g,%u\n",l,mlmc_p.Nl[l],Nl_scale[l],time,sim_counter);
    }
    for (l=0;l<=mlmc_p.L;l++)
    {
        fprintf(stderr,"NACC = %d\n",abc_p.nacc);
        mlmc_p.Nl[l] = Nl_scale[l]*abc_p.nacc;
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
    fprintf(stdout,"\"theta\",\"F\",\"SEC\",\"NSIMS\"\n");
    for (i=0;i<M.G.numPoints;i++)
    {
        fprintf(stdout,"%g,%g,%g,%u\n",M.G.coords[i],M.F[i],time,sim_counter);
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



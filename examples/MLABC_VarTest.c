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

/* Example of ABC-PCR samples for TB transmission model from Tanaka et al.*/
ABC_Parameters abc_p;
MLMC_Parameters mlmc_p;


Dataset data; /*not really used except as a template for simulated data*/
/*generate grid for indicator thresholds*/
ndgrid GenNdGrid(int,double *,double *,int,double *);
/*discrepency metric */
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
    if (support == NULL)
    {
        theta[0] = durngus(0.0,5.0);/*alpha: birth rate*/
        theta[1] = durngus(0.0,theta[0]);/*delta: death rate*/
        theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
        while (theta[2] < 0.0)
        {
            theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
        }
    }
    else /*for now we only control alpha and delta support*/
    {
        theta[0] = durngus(support[0],support[3]);
        theta[1] = (support[4] < theta[0]) ? durngus(support[1],support[4]) : durngus(support[1],theta[0]);
        theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
        while (theta[2] < support[2] || theta[2] > support[5])
        {
            theta[2] = durngns(0.198,0.06735);/*theta: muation rate*/
        }
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

double deltas[3];

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
    double *gH_d;
    CDF_estimate M;
    clock_t start_t, end_t;
    double time;
    unsigned int i,l,t;
    unsigned int numTrials;
    double *Yl,*Ylu,*Fl,*EYl,*EYlu,*EFl,*EYl2,*EYlu2,*EFl2;
    /*only 2 args*/
    if (argc  < 4)
    {
        fprintf(stderr,"Usage: %s N eps numTrials\n",argv[0]);
    }
    else
    {
        /*set up ABC params*/
        mlmc_p.presample_trials = (unsigned int)atoi(argv[1]);
        abc_p.eps = (double)atof(argv[2]);
        abc_p.nmax = 0;
        abc_p.k = 3;
        abc_p.support = (double*)malloc(6*sizeof(double));
        abc_p.support[0] = 0.0;
        abc_p.support[1] = 0.0;
        abc_p.support[2] = 0.0;
        abc_p.support[3] = 5.0;
        abc_p.support[4] = 5.0;
        abc_p.support[5] = 0.5;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
        numTrials = (unsigned int)atoi(argv[3]);
    }
 
    /*set up MLMC parameters*/
    mlmc_p.L = 10 ;
    mlmc_p.eps0 = 1.0;
    mlmc_p.eps_l = (double*)malloc((mlmc_p.L+1)*sizeof(double));
    mlmc_p.Nl = (unsigned int*)malloc((mlmc_p.L+1)*sizeof(unsigned int));
    Nl_scale = (double*)malloc((mlmc_p.L+1)*sizeof(double));
    mlmc_p.eps_l[0] = 1.0;
    mlmc_p.eps_l[mlmc_p.L] = abc_p.eps;
    for (i=1;i<mlmc_p.L;i++)
    {
        mlmc_p.eps_l[i] = mlmc_p.eps_l[mlmc_p.L] + 0.5*(mlmc_p.eps_l[i-1] - mlmc_p.eps_l[mlmc_p.L]);
    }
    mlmc_p.target_RMSE = 0.1;
    mlmc_p.presample = 1; 

    /*initialise grid*/
    deltas[0] = (abc_p.support[3] - abc_p.support[0])/99.0; 
    deltas[1] = (abc_p.support[4] - abc_p.support[1])/99.0; 
    deltas[2] = (abc_p.support[5] - abc_p.support[2])/99.0; 
    M.G = GenNdGrid(3,abc_p.support,abc_p.support+3,100,deltas);
    M.F = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.V = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.g = &g; 

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



    Yl = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    Ylu = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    Fl = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    
    EYl = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    EYlu = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    EFl = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    EYl2 = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    EYlu2 = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    EFl2 = (double*)malloc((mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    
    memset(EYl,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    memset(EYlu,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    memset(EFl,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    memset(EYl2,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    memset(EYlu2,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
    memset(EFl2,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));

    for (t=0;t<numTrials;t++)
    {
        memset(Yl,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
        memset(Ylu,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
        memset(Fl,0,(mlmc_p.L+1)*(M.G.numPoints)*sizeof(double));
        start_t = clock();
        //dmlabcnls(abc_p,mlmc_p.presample_trials,mlmc_p,&data,&M,mlmc_p.Nl,Nl_scale);
        dmlabcnlps(abc_p,mlmc_p.presample_trials,mlmc_p,&data,&M,Yl,Ylu,Fl);
        for (i=0;i<(mlmc_p.L+1)*M.G.numPoints;i++){
            EYl[i] += Yl[i];
            EYlu[i] += Ylu[i];
            EFl[i] += Fl[i];
            EYl2[i] += Yl[i]*Yl[i];
            EYlu2[i] += Ylu[i]*Ylu[i];
            EFl2[i] += Fl[i]*Fl[i];
        }
        end_t = clock();

        time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
        fprintf(stderr,"\nSim %d or %d time = %g secs.\n",t+1,numTrials,time);
    }
    for (i=0;i<(mlmc_p.L+1)*M.G.numPoints;i++){
        EYl[i] /= (double)numTrials;
        EYlu[i] /= (double)numTrials;
        EFl[i] /= (double)numTrials;
        EYl2[i] /= (double)numTrials;
        EYlu2[i] /= (double)numTrials;
        EFl2[i] /= (double)numTrials;
    }
    for (i=0;i<(mlmc_p.L+1)*M.G.numPoints;i++){
        EYl2[i] -= EYl[i]*EYl[i];
        EYlu2[i] -= EYlu[i]*EYlu[i];
        EFl2[i] -= EFl[i]*EFl[i];
    }

    /*write results*/
    fprintf(stdout,"\"Level\",\"V Y_l\",\"V Y_l uc\",\"V F_l sl\"\n");
    for (l=0;l<mlmc_p.L+1;l++){
        double VYl,VYlu,VFl;
        VYl = EYl2[l*M.G.numPoints];
        for (i=0;i<M.G.numPoints;i++){
            VYl = (VYl < EYl2[l*M.G.numPoints + i]) ?  EYl2[l*M.G.numPoints + i] : VYl;
        }
        VYlu = EYlu2[l*M.G.numPoints];
        for (i=0;i<M.G.numPoints;i++){
            VYlu = (VYlu < EYlu2[l*M.G.numPoints + i]) ?  EYlu2[l*M.G.numPoints + i] : VYlu;
        }
        VFl = EFl2[l*M.G.numPoints];
        for (i=0;i<M.G.numPoints;i++){
            VFl = (VFl < EFl2[l*M.G.numPoints + i]) ?  EFl2[l*M.G.numPoints + i] : VFl;
        }
        fprintf(stdout,"%d, %g,%g,%g\n",l,VYl,VYlu,VFl);
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



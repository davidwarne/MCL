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
            theta[i*3] = durngus(0.0,28.0); 
            theta[i*3 + 1] = durngus(0.0,0.04); 
            theta[i*3 + 2] = durngus(0.0,28.0); 
        }
    }
    else
    {
        for (i=0;i<nsamples;i++)
        {
            theta[i*3] = durngus(support[0],support[3]); 
            theta[i*3 + 1] = durngus(support[1],support[4]); 
            theta[i*3 + 2] = durngus(support[2],support[5]); 
        }
    }
    return 0;
}

#define NUMREAL 3
/*simulation*/
int simulate(void *sim,double * theta, Dataset * data_s)
{
    unsigned int n,m,nt,i,j;
    double *Y_r;
    double *X_r;
    double Y0[3] = {1,1000,1000};
    double nu[9] = {0,1,0,0,-1,1,0,0,-1};
    double nu_minus[9] = {1,1,0,0,1,1,0,0,1}; 
    double T[19] = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0}; 
    unsigned int d[2] = {1,2};

    /*intial conditions*/
    Y_r = (double*)data_s->fields[0].data_array;
    nt = data_s->fields[0].numRows;
    n = data_s->fields[0].numCols;

    for (j=0;j<n*nt;j++)
    {
        Y_r[j] = 0;
    }

    if (n != 2 || nt != 19)
    {
        fprintf(stderr,"Fatal Error: Invalid dataset [n,nt = %d %d].\n",n,nt);
        exit(1);
    }

    m = 3; /*params are c1,c2,c3*/
    
    for (i=0;i<NUMREAL;i++)
    {
        X_r = (double*)data_s->fields[i+1].data_array;
        degils(m,n+1,nt,T,Y0,nu_minus,nu,theta,n,d,X_r);
        for (j=0;j<n*nt;j++)
        {
            Y_r[j] += X_r[j];
        }
    }
    
    for (j=0;j<n*nt;j++)
    {
        Y_r[j] /= ((double)NUMREAL);
    }

    sim_counter++;
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
        abc_p.eps = 1800.0;
        abc_p.nmax = 0;
        abc_p.k = 3;
        abc_p.support = (double*)malloc(6*sizeof(double));
        abc_p.support[0] = 0.0;
        abc_p.support[1] = 28.0;
        abc_p.support[2] = 0.0;
        abc_p.support[3] = 0.04;
        abc_p.support[4] = 0.0;
        abc_p.support[5] = 28.0;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
 
    /*set up SMC parameters*/
    mlmc_p.L = 4;
    mlmc_p.eps_l = (double*)malloc(5*sizeof(double));
    mlmc_p.eps_l[0] = 4000.0;
    mlmc_p.eps_l[1] = 2900.0; 
    mlmc_p.eps_l[2] = 2000.0;
    mlmc_p.eps_l[3] = 1900.0;
    mlmc_p.eps_l[4] = 1800.0;
    mlmc_p.Nl = (unsigned int*)malloc((mlmc_p.L+1)*sizeof(unsigned int));
    Nl_scale = (double*)malloc((mlmc_p.L+1)*sizeof(double));
    mlmc_p.presample = 1; 
    mlmc_p.presample_trials = 100; 
    
    /*initialise grid: we define it only in */
    deltas[0] = (abc_p.support[3] - abc_p.support[0])/99.0; 
    deltas[1] = (abc_p.support[4] - abc_p.support[1])/99.0; 
    deltas[2] = (abc_p.support[5] - abc_p.support[2])/99.0; 
    M.G = GenNdGrid(abc_p.k,abc_p.support,abc_p.support+abc_p.k,100,deltas);
    M.F = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.V = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.g = &g; 

    /*initialise our RNG library*/
    SSAL_Initialise(argc,argv);

    /*build data template data + space*/
    data.numFields = 1 + NUMREAL ;
    data.fields = (field *)malloc(data.numFields*sizeof(field));
    for (i=0;i<data.numFields;i++)
    {
        data.fields[i].numCols = 2;
        data.fields[i].numRows = 19;
        data.fields[i].type = REAL64_DATA;
        data.fields[i].data_array = malloc(38*sizeof(double));
        data.fields[i].numBytes = 38*sizeof(double);
    }
    /* populate data (pre-computed)*/
    Y_d = (double*)data.fields[0].data_array;
  Y_d[0]= 1026.6; Y_d[19]= 1084.7;
  Y_d[1]=  1007.9;Y_d[20]= 1053.9;
  Y_d[2]=  998.5; Y_d[21]= 1050.6;
  Y_d[3]=  876.6; Y_d[22]= 894.4;
  Y_d[4]=  1064.7;Y_d[23]= 1258.5;
  Y_d[5]=  815.7; Y_d[24]= 640.8;
  Y_d[6]=  1047.0;Y_d[25]= 1686.8;
  Y_d[7]=  780.0; Y_d[26]= 529.8;
  Y_d[8]=  1156.1;Y_d[27]= 1937.6;
  Y_d[9]=  631.1; Y_d[28]= 449.0;
  Y_d[10]=  1724.0;Y_d[29]= 1177.3;
  Y_d[11]=  564.0; Y_d[30]= 560.1;
  Y_d[12]=  1424.8;Y_d[31]= 1371.0;
  Y_d[13]=  377.2; Y_d[32]= 815.5;
  Y_d[14]=  2028.0;Y_d[33]= 637.6;
  Y_d[15]=  812.1; Y_d[34]= 1338.8;
  Y_d[16]=  903.0; Y_d[35]= 287.2;
  Y_d[17]=  1575.3;Y_d[36]= 1631.3;
  Y_d[18]=  466.8; Y_d[37]= 914.6;
    
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
    fprintf(stdout,"\"c1\",\"c2\",\"c3\",\"F\",\"SEC\",\"NSIMS\"\n");
    for (i=0;i<M.G.numPoints;i++)
    {
        fprintf(stdout,"%g,%g,%g,%g,%g,%u\n",M.G.coords[i*3],M.G.coords[i*3+1],M.G.coords[i*3+2],M.F[i],time,sim_counter);
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



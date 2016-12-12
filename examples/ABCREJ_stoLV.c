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
        theta[i*3] = durngus(0.0,28.0); 
        theta[i*3 + 1] = durngus(0.0,0.04); 
        theta[i*3 + 2] = durngus(0.0,28.0); 
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
        abc_p.eps = 1800.0;
        abc_p.nmax = 0;
        abc_p.k = 3;
        abc_p.support = NULL;
        abc_p.sim = NULL;
        abc_p.rho = &rho;
        abc_p.p = &prior;
        abc_p.pd = NULL;
        abc_p.s = &simulate;
    }
 
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
    theta = (double *)malloc(abc_p.nacc*abc_p.k*sizeof(double));

    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC-SMC*/
    start_t = clock();
    dabcrs(abc_p,&data,theta,NULL);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stdout,"\"Sample\",\"c1\",\"c2\",\"c3\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%g,%u,%u\n",i,theta[i*3],theta[i*3+1],theta[i*3+2],time,abc_p.nacc,sim_counter);
    }

    return 0;
}


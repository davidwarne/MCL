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
/* Example of ABC-rejection samples for a deteministic repressilator system from Toni et al.*/
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
    return sqrt(d);
}

/*prior sampler a ~ U(-10,10) b ~ U(-10,10)*/
int prior(unsigned int dim, unsigned int nsamples, double* support,double * theta)
{
    unsigned int i,j;
    for (i=0;i<nsamples;i++)
    {
        theta[i*4] = durngus(-2.0,10.0); /*\alpha_0*/
        theta[i*4+1] = durngus(0,10.0); /*n*/
        theta[i*4+2] = durngus(-5.0,20.0); /*\beta*/
        theta[i*4+3] = durngus(500.0,2500.0); /*\alpha*/
    }
    return 0;
}

/* Deterministic Repressilator Genetic regulatory network
 *
 */
void detRep(double *Y, unsigned int d, double *params,unsigned int m, double t,double *f_r)
{
   double alpha_0,n,beta,alpha;
   double m_1,p_1,m_2,p_2,m_3,p_3;
   alpha_0 = params[0];
   n = params[1];
   beta = params[2];
   alpha = params[3];
   /* m_i mRNA and p_i produced protien*/
   m_1 = Y[0];
   p_1 = Y[1];
   m_2 = Y[2];
   p_2 = Y[3];
   m_3 = Y[4];
   p_3 = Y[5];

   f_r[0] = -m_1 + (alpha)/(1.0 + pow(p_3,n)) + alpha_0;
   f_r[1] = -beta*(p_1 - m_1);
   f_r[2] = -m_2 + (alpha)/(1.0 + pow(p_1,n)) + alpha_0;
   f_r[3] = -beta*(p_2 - m_2);
   f_r[4] = -m_3 + (alpha)/(1.0 + pow(p_2,n)) + alpha_0;
   f_r[5] = -beta*(p_3 - m_3);
}

#define TSTEP 5e-3
/*simulation*/
int simulate(void *sim,double * theta, Dataset * data_s)
{
    unsigned int n,m,nt,i;
    double *Y_r;
    double Y0[6] = {0.0,2.0,0.0,1.0,0.0,3.0};
    double T[12] = {0.6,4.2,6.2,8.6,13.4,16.0,21.4,27.6,34.4,39.8,40.6,45.2}; 
    unsigned int d[3] = {0,2,4};
    double h;
    h = TSTEP;

    /*intial conditions*/
    Y_r = (double*)data_s->fields[0].data_array;
    nt = data_s->fields[0].numRows;
    n = data_s->fields[0].numCols;

    if (n != 3 || nt != 12)
    {
        fprintf(stderr,"Fatal Error: Invalid dataset [n,nt = %d %d].\n",n,nt);
        exit(1);
    }

    m = 4; /*params are  alpha_0,n,beta and alpha*/
    
    /*solve ODE with RK4 method*/
    drk4s(m,n*2,nt,T,theta,Y0,&detRep,n,d,h,Y_r);

    sim_counter++;
    return 0;
}

int main(int argc, char ** argv)
{
    double *theta,*r;
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
        abc_p.eps = 35.0;
        abc_p.nmax = 0;
        abc_p.k = 4;
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
    data.fields[0].numRows = 12;
    data.fields[0].type = REAL64_DATA;
    data.fields[0].data_array = malloc(36*sizeof(double));
    data.fields[0].numBytes = 36*sizeof(double);
   
    /* populate data (pre-computed)*/
    Y_d = (double*)data.fields[0].data_array;
    Y_d[0]= 2.04;   Y_d[12]=  28.99; Y_d[24]= 20.96;
    Y_d[1]= 32.19;  Y_d[13]= 11.29;  Y_d[25]= 7.49;
    Y_d[2]= 4.13;   Y_d[14]= 10.61;  Y_d[26]= 44.25;
    Y_d[3]= 2.15;   Y_d[15]= 55.27;  Y_d[27]= 7.12;
    Y_d[4]= 5.09;   Y_d[16]= 9.49;   Y_d[28]= 60.52;
    Y_d[5]= 1.07;   Y_d[17]= 68.56;  Y_d[29]= 8.10;
    Y_d[6]= 3.67;   Y_d[18]= 10.62;  Y_d[30]= 63.76;
    Y_d[7]= 39.01;  Y_d[19]= -1.95;  Y_d[31]= 22.90;
    Y_d[8]= 73.83;  Y_d[20]= 3.53;   Y_d[32]= 6.27;
    Y_d[9]= 8.54;   Y_d[21]= 63.87;  Y_d[33]= 10.59;
    Y_d[10]= 17.62;  Y_d[22]= 39.68;  Y_d[34]= 6.50;
    Y_d[11]= 11.96;  Y_d[23]= -0.60;  Y_d[35]= 70.56;
    /*allocate output array*/
    theta = (double *)malloc(abc_p.nacc*abc_p.k*sizeof(double));
    r= (double *)malloc(abc_p.nacc*sizeof(double));

    /*test smallest eps*/
    //{
    //    Dataset *data_s;
    //    data_s = copyDataset(&data); 
    //    theta[0] = 1.0;
    //    theta[1] = 2.0;
    //    theta[2] = 5.0;
    //    theta[3] = 1000.0;
    //    simulate(NULL,theta,data_s);
    //    fprintf(stderr,"eps = %f\n",rho(&data,data_s));
    //}
    sim_counter = 0; /*for performance metric*/
    /*generate posterior samples using ABC-MCMC*/
    start_t = clock();
    dabcrs(abc_p,&data,theta,r);
    end_t = clock();
    time = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);
    /*write output*/
    fprintf(stdout,"\"Sample\",\"alpha_0\",\"n\",\"beta\",\"alpha\",\"SEC\",\"NACC\",\"NSIMS\"\n");
    for (i=0;i<abc_p.nacc;i++)
    {
        fprintf(stdout,"%u,%g,%g,%g,%g,%g,%u,%u\n",i,theta[i*4],theta[i*4+1],theta[i*4+2],theta[i*4+3],time,abc_p.nacc,sim_counter);
    }

    return 0;
}


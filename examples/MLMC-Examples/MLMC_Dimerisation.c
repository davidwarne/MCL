/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2017  David J. Warne
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

/**
 * @file MLMC_Dimerisation.c 
 *
 * @brief An examples of Biased and Exact MLMC 
 * @details Application is a simple 2nd order chemical reaction network 
 * modelling protien dimerisation as per Lester et al[1].
 * 
 * @note [1] C Lester, RE Baker, MB Giles and CA Yates (2016) Extending the 
 * Multi-level method fo the simulation of stochastic biological systems. 
 * Bulletin of Mathematica Biology. 78:1640--1677.
 */

/* values of inverse normal CDF for confidence intervals*/
#define CI_95 1.96
#define CI_99 2.5758

#define CI CI_99

SSAL_CRN dimer_model;
int EXACT;
double
f(double X)
{
    double sigma, Y;
    sigma = 100;
    Y = 3708;
    return 1e7*exp(-(Y-X)*(Y-X)/(2.0*sigma*sigma))/(sigma*sqrt(2.0*M_PI));
}


int 
simulate_sl(int n,double * params, double *D)
{
    double T;
    double p;
    T = params[1];
    p = params[3];
    for (int i = 0;i<n;i++)
    {
       // degils(dimer_model.M,dimer_model.N,1,&T,dimer_model.X0,
       //        dimer_model.nu_minus,dimer_model.nu,dimer_model.c,
       //        dimer_model.nvar,dimer_model.vars,D+i);
       demnrms(dimer_model.M,dimer_model.N,1,&T,dimer_model.X0,
              dimer_model.nu_minus,dimer_model.nu,dimer_model.c,
              dimer_model.nvar,dimer_model.vars,D+i);
            //D[i] = pow(D[i],p);
            D[i] = f(D[i]);
    }
}


int
simulate_ml(int l, int nl, int K, double *params, double *D)
{
    double tau0;
    int M;
    int E;
    double T;
    double p;
    /*extract parameters*/
    tau0 = params[0];
    T = params[1];
    M = (int)params[2];
    p = params[3];
    if (K != 1 || dimer_model.nvar != 1)
    {
        return 1;
    }

    if (l == 0) /* coarse level*/
    {
        for (int i=0;i<nl;i++)
        {
            datauls(dimer_model.M,dimer_model.N,1,&T,dimer_model.X0,
                    dimer_model.nu_minus,dimer_model.nu,dimer_model.c,
                    dimer_model.nvar,dimer_model.vars,tau0, D + i);
            //D[i] = pow(D[i],p);
            D[i] = f(D[i]);
        }
    }
    else if (l == EXACT)
    {
        double taul = tau0/pow((double)M,(double)(l-1));
        double D_e, D_c;
        for (int i=0;i<nl;i++)
        {
            daectauls(dimer_model.M,dimer_model.N,1,&T,dimer_model.X0,
                    dimer_model.nu_minus,dimer_model.nu,dimer_model.c,
                    dimer_model.nvar,dimer_model.vars,taul, &D_e,&D_c);
            //D[i] = pow(D_e,p) - pow(D_c,p);
            D[i] = f(D_e) - f(D_c);
        }

    }
    else /* coupled tau-leap*/
    {
        double taul = tau0/pow((double)M,(double)l);
        double D_f, D_c;
        for (int i=0;i<nl;i++)
        {
            dactauls(dimer_model.M,dimer_model.N,1,&T,dimer_model.X0,
                    dimer_model.nu_minus,dimer_model.nu,dimer_model.c,
                    dimer_model.nvar,dimer_model.vars,taul,M, &D_f,&D_c);
            //D[i] = pow(D_f,p) - pow(D_c,p);
            D[i] = f(D_f) - f(D_c);
        }

    }
    
    return 0;
}

int
main(int argc, char ** argv)
{
    /* Target standard deviation of estimator (i.e. the standard error) */
    double eps;     
    double mol; /*Target confinence interval*/
    double params[5];
    int L;
    int *nl;
    int n;
    double *D;
    double ED;
    double sig2;
    unsigned int ml;
    /* Init the RNG*/
    SSAL_Initialise(argc,argv);
    
    /* import Dimerisation model*/
    dimer_model = SSAL_ImportLSBML("dimer.json");
    /*Sanity check*/ 
    SSAL_WriteChemicalReactionNetwork(stderr,dimer_model);
    
    /*we only want to track the dimers*/
    dimer_model.nvar = 1;
    dimer_model.vars[0] = 2; /* we have [M,P,D]*/

    ml = 1;
    if (argc < 6)
    {
        fprintf(stderr,"Usage: %s moment mol tau0 L M \n",argv[0]);
        params[3] = (double)(atoi(argv[1]));
        mol = (double)(atof(argv[2]));
        ml = 0;
    }
    else
    {
        params[0] = (double)atof(argv[3]);
        L = (int)atoi(argv[4]);
        params[2] = (double)(atoi(argv[5]));
        params[3] = (double)(atoi(argv[1]));
        mol = (double)(atof(argv[2]));
    }
    eps = mol/CI;

    params[1] = 1.0; 
    if (ml == 1)
    {
        nl = (int *)malloc((L+1)*sizeof(int));
        EXACT = L;
        /*estimate nl's*/
        dmlmcnls(L,100,eps,params,&simulate_ml,dimer_model.nvar,nl);  
        for (int l=0;l<=L;l++)
        {
            fprintf(stderr,"n_%d = %d\n",l,nl[l]);
        }
        /*run MLMC*/
        dmlmcs(L,nl,params,&simulate_ml,dimer_model.nvar,&ED,&sig2);

        fprintf(stdout,"\"tau0\",\"L\",\"M\",\"eps\",\"E[D]\",\"sig\"\n");
        fprintf(stdout,"%g,%d,%d,%g,%g+-%g,%g\n",params[0],L,(int)params[2],eps,
                ED,CI*sqrt(sig2),sig2);
    }
    else
    {
        /*Compare with Standard Monte Carlo*/
        n = 100;
        D = malloc(dimer_model.N*n*sizeof(double));
        simulate_sl(n,params,D);
        dmcintv(n,1,D,&ED,&sig2);
        n = (int)ceil(sig2/(eps*eps));
        fprintf(stderr,"n = %d\n",n);
        D = (double*)realloc(D, dimer_model.N*n*sizeof(double));
        simulate_sl(n,params,D);
        dmcintv(n,1,D,&ED,&sig2);
        sig2 /= (double)n;
        fprintf(stdout,"\"eps\",\"E[D]\",\"sig\"\n");
        fprintf(stdout,"%g,%g+-%g,%g\n",eps,ED,CI*sqrt(sig2),sig2);

    }
}

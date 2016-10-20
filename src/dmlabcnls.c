/* LIBABC : A Computational Library for approximate Bayesian computation.
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
/**
 * @brief samples runtime and variance for each bias correction term in ML estimater to
 * select optimal N_l based on these.
 */
int dmlabcnls(ABC_Parameters abc_p, int trials,ML_Parameters ml_p, Dataset *data, CDF_estimate *M, unsigned int *Nl)
{
    int l,j,i,k;
    ABC_Parameters *abc_pl;
    double *E,*V;
    double *times;
    double *sigma2;
    double *thetal, *thetalm1;
    double *rhol,*rholm1;
    double c;

    if (ml_p.presample == 0)
    {
        double gamma,beta,alpha,K,J,h;
        gamma = ml_p.gamma;
        beta = ml_p.beta;
        alpha = ml_p.alpha;
        K = (double)(ml_p.K);
        J = (double)(M->G.numPoints);
        h = ml_p.target_RMSE; 

        if (gamma < beta)
        {
        }
        else if (gamma == beta)
        {
        }
        else
        {
            
        }

        return 0;
    }

    /* allocate sequence of ABC parameters*/
    abc_pl = (ABC_Parameters*)malloc((ml_p.L+1)*sizeof(ABC_Parameters));
    for (l=0;l<=ml_p.L;l++)
    {
        abc_pl[l] = abc_p;
        abc_pl[l].eps = ml_p.eps_l[l];
        abc_pl[l].nacc = trials;
        abc_pl[l].support = (SSAL_real_t *)malloc(abc_p.k*2*sizeof(double));
        memcpy(abc_pl[l].support,abc_p.support,abc_p.k*2*sizeof(SSAL_real_t));

    }
    
    E = (double *)malloc(M->G.numPoints*sizeof(double));
    V = (double *)malloc(M->G.numPoints*sizeof(double));

    sigma2 = (double*)malloc((ml_p.L+1)*sizeof(double));
    times = (double*)malloc((ml_p.L+1)*sizeof(double));

    thetal = (double*)malloc(trials*abc_p.k*sizeof(double));
    thetalm1 = (double*)malloc(trials*abc_p.k*sizeof(double));
    rhol = (double*)malloc(trials*abc_p.k*sizeof(double));

    /*build temporary */

    for (j=0;j<M->G.numPoints;j++)
    {
        M->F[j] = 0 ;
    }
    /*begin timing trials*/
    for (l=0;l<=ml_p.L;l++)
    {
        clock_t start_t,end_t;   

        int temp_nacc;
        double temp_support[255];

        start_t = clock();
        if (l == 0)
        {
            /* sample ABC  inital biased estimate*/
            dabcrs(abc_pl[l], data,thetal, rhol);
            /* Integrate smoothed indicator function */
            for (j=0;j<M->G.numPoints;j++)
            {
                dmcint(abc_pl[l].nacc,abc_pl[l].k,thetal,M->G.coords + j*M->G.dim,M->g,E+j,V+j);
            }

        }
        else
        {
            memcpy(thetalm1,thetal,abc_pl[l-1].k*abc_pl[l-1].nacc*sizeof(double));
            memcpy(abc_pl[l].support,thetalm1,abc_pl[l].k*sizeof(double));
            memcpy(abc_pl[l].support+abc_pl[l].k,thetalm1,abc_pl[l].k*sizeof(double));
            /*update support for sampling bias term */
            for (i=1;i<abc_pl[l-1].nacc;i++)
            {
                for (j=0;j<abc_pl[l-1].k;j++)
                {
                    if (thetalm1[i*abc_pl[l].k +j] < abc_pl[l].support[j])
                    {
                        abc_pl[l].support[j] = thetalm1[i*abc_pl[l].k + j];
                    }

                    if (thetalm1[i*abc_pl[l].k +j] > abc_pl[l].support[abc_pl[l].k + j])
                    {
                        abc_pl[l].support[abc_pl[l].k + j] = thetalm1[i*abc_pl[l].k + j];
                    }
                }
            }
            /*sample coupled levels for bias correction terms*/
            dabcrjs(abc_pl[l-l],abc_pl[l], data, M,thetalm1,thetal,rhol);
            /* integrate the difference in smoothed indicator functions*/
            for (j=0;j<M->G.numPoints;j++)
            {
               dmcintd(abc_pl[l].nacc,abc_pl[l].k,thetal,thetalm1,M->G.coords +j*M->G.dim,M->g,E+j,V+j);
            }
        }
        end_t = clock();

        /* compute Var[g_sj] = |(E[g_sj^2] - E[g_sj])|_\infty^2*/
        sigma2[l] = V[0];
        for (j=0;j<M->G.numPoints;j++)
        {
            if (V[j] > sigma2[l])
            {
                sigma2[l] = V[j];
            }
        }
        times[l] = ((double)(end_t - start_t))/((double)CLOCKS_PER_SEC);

        for (j=0;j<M->G.numPoints;j++)
        {
            M->F[j] += E[j];
        }
    }

    /*use optimal solution based onn Lagrange multiplier*/
    c = 0;
    for (l=0;l<=ml_p.L;l++)
    {
        c += sqrt(sigma2[l]*times[l]);
    }
    for (l=0;l<=ml_p.L;l++)
    {
        Nl[l] = (unsigned int)ceil(log((double)M->G.numPoints)*c*(sqrt(sigma2[l])/sqrt(times[l]))/(ml_p.target_RMSE*ml_p.target_RMSE));
    }

    for (l=0;l<=ml_p.L;l++)
    {
        Nl[l] = (Nl[l] < trials) ? trials : Nl[l];
        printf("N[%d] = %d, Var[%d] = %g, eps[%d] = %g,c[%d] = %g\n",l,Nl[l],l,sigma2[l],l,ml_p.eps_l[l],l,times[l]);
    }

    free(times);
    free(sigma2);
    free(thetal);
    free(thetalm1);
    free(rhol);
    free(abc_pl);
    return 0;
}

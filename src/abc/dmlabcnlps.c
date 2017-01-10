/* MCL: Monte Carlo Library.
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
#include <time.h>

#define MINCLOCKS 5000
/**
 * @brief runs a fixed set of sampled per level and returns each term in both correlated and uncorrelated versions AND as a single level version.
 */
int dmlabcnlps(ABC_Parameters abc_p, int trials,MLMC_Parameters ml_p, Dataset *data, CDF_estimate *M, double *Yl, double *Ylu, double *Fl)
{
    int l,j,i,k;
    ABC_Parameters *abc_pl;
    double *E,*V;
    double *thetal, *thetalm1;
    double *rhol,*rholm1;
    double c;
    CDF_estimate *M_mon;
    M_mon = (CDF_estimate*)malloc(sizeof(CDF_estimate));
    M_mon->G = M->G;
    M_mon->F = (double*)malloc(M_mon->G.numPoints*sizeof(double));
    M_mon->V = (double*)malloc(M_mon->G.numPoints*sizeof(double));
    M_mon->g = M->g;
    
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

    thetal = (double*)malloc(trials*abc_p.k*sizeof(double));
    thetalm1 = (double*)malloc(trials*abc_p.k*sizeof(double));
    rhol = (double*)malloc(trials*abc_p.k*sizeof(double));

    /*build temporary */

    for (j=0;j<M->G.numPoints;j++)
    {
        M->F[j] = 0 ;
    }
    
    fprintf(stderr,"\nRunning correlated terms\n");
    /*begin timing trials*/
    for (l=0;l<=ml_p.L;l++)
    {

        int temp_nacc;
        double temp_support[255];
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
            dabcrjs(abc_pl[l-l],abc_pl[l], data, M_mon,thetalm1,thetal,rhol);
            /* integrate the difference in smoothed indicator functions*/
            for (j=0;j<M->G.numPoints;j++)
            {
               dmcintd(abc_pl[l].nacc,abc_pl[l].k,thetal,thetalm1,M->G.coords +j*M->G.dim,M->g,E+j,V+j);
            }
        }
        
        for (j=0;j<M->G.numPoints;j++)
        {
            M->F[j] += E[j];
        }
        monotone(M,M_mon);
        
        for (j=0;j<M->G.numPoints;j++)
        {
            Yl[l*(M->G.numPoints) +j] = E[j];
        }
        fprintf(stderr," %d",l);
    }
    /*reset for uncorrelated trial*/
    for (l=0;l<=ml_p.L;l++)
    {
        abc_pl[l] = abc_p;
        abc_pl[l].eps = ml_p.eps_l[l];
        abc_pl[l].nacc = trials;
        abc_pl[l].support = (SSAL_real_t *)malloc(abc_p.k*2*sizeof(double));
        memcpy(abc_pl[l].support,abc_p.support,abc_p.k*2*sizeof(SSAL_real_t));

    }
    
    for (j=0;j<M->G.numPoints;j++)
    {
        M->F[j] = 0 ;
    }
    
    fprintf(stderr,"\nRunning uncorrelated terms\n");
    /*begin timing trials*/
    for (l=0;l<=ml_p.L;l++)
    {

        int temp_nacc;
        double temp_support[255];
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
            dabcrs(abc_pl[l], data,thetal, rhol);
            dabcrs(abc_pl[l-1], data,thetalm1, NULL);
            /* integrate the difference in smoothed indicator functions*/
            for (j=0;j<M->G.numPoints;j++)
            {
               dmcintd(abc_pl[l].nacc,abc_pl[l].k,thetal,thetalm1,M->G.coords +j*M->G.dim,M->g,E+j,V+j);
            }
        }
        
        for (j=0;j<M->G.numPoints;j++)
        {
            Ylu[l*(M->G.numPoints) +j] = E[j];
        }
        fprintf(stderr," %d",l);
    }
    /*reset for single level trial*/
    for (l=0;l<=ml_p.L;l++)
    {
        abc_pl[l] = abc_p;
        abc_pl[l].eps = ml_p.eps_l[l];
        abc_pl[l].nacc = trials;
        abc_pl[l].support = (SSAL_real_t *)malloc(abc_p.k*2*sizeof(double));
        memcpy(abc_pl[l].support,abc_p.support,abc_p.k*2*sizeof(SSAL_real_t));

    }
    
    fprintf(stderr,"\nRunning single level terms\n");
    for (l=0;l<=ml_p.L;l++)
    {

        int temp_nacc;
        double temp_support[255];
        /* sample ABC  inital biased estimate*/
        dabcrs(abc_pl[l], data,thetal, rhol);
        /* Integrate smoothed indicator function */
        for (j=0;j<M->G.numPoints;j++)
        {
            dmcint(abc_pl[l].nacc,abc_pl[l].k,thetal,M->G.coords + j*M->G.dim,M->g,E+j,V+j);
        }
        
        for (j=0;j<M->G.numPoints;j++)
        {
            Fl[l*(M->G.numPoints) +j] = E[j];
        }
        fprintf(stderr," %d",l);
    }

    free(thetal);
    free(thetalm1);
    free(rhol);
    free(abc_pl);
    return 0;
}


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

/**
 * @brief Generate esitmates for the the Multivariate distribution function
 * at given grid points using a multilevel ABC method.
 *
 * @param sim and SSAL simulation struct
 * @param abc_p ABC Parameters
 * @param ml_p Multilevel Parameters
 * @param data dataset to perform inference on
 * @param M the resulting estimated CDF 
 *
 */
int 
dmlabccdfs(ABC_Parameters abc_p,MLMC_Parameters ml_p, Dataset *data, 
           CDF_estimate *M)
{
    /*initialise*/
    double *E;
    double *V;
    ABC_Parameters *abc_pl;
    size_t max_n;
    double *thetal, *thetalm1;
    double *rhol;
    unsigned int i,j,l,k;
   
    CDF_estimate *M_mon;
    M_mon = (CDF_estimate*)malloc(sizeof(CDF_estimate));
    M_mon->G = M->G;
    M_mon->F = (double*)malloc(M_mon->G.numPoints*sizeof(double));
    M_mon->V = (double*)malloc(M_mon->G.numPoints*sizeof(double));
    M_mon->g = M->g;
    max_n = 0;
    for (l=0;l<=ml_p.L;l++)
    {
        max_n = (ml_p.Nl[l] > max_n) ? ml_p.Nl[l] : max_n;
    }
    /*parameters for each level*/
    if ((abc_pl = (ABC_Parameters *)malloc((ml_p.L+1)*sizeof(ABC_Parameters))) 
         ==  NULL)
    {
        return 1;
    }
    for (l=0;l<=ml_p.L;l++)
    {
        abc_pl[l] = abc_p;
        abc_pl[l].eps = ml_p.eps_l[l];
        abc_pl[l].support = (SSAL_real_t *)malloc(abc_p.k*2*sizeof(double));
        memcpy(abc_pl[l].support,abc_p.support,abc_p.k*2*sizeof(SSAL_real_t));
        abc_pl[l].nacc = ml_p.Nl[l];
    }

    l=0;
    /*allocate expectations*/
    if ((E = (double*) malloc(M->G.numPoints*sizeof(double))) == NULL)
    {
        return 1;
    }
    if ((V = (double*) malloc(M->G.numPoints*sizeof(double))) == NULL)
    {
        return 1;
    }

    if ((thetal = (double *)malloc(((size_t)max_n)*((size_t)abc_pl[0].k)*
                                   sizeof(double))) == NULL)
    {
        
        return 1;
    }
    if ((thetalm1 = (double *)malloc(max_n*abc_pl[0].k*sizeof(double))) == NULL)
    {
        return 1;
    }
    if ((rhol = (double *)malloc(max_n*sizeof(double))) == NULL)
    {
        return 1;
    }
   
 
        /*generate N_0 samples*/
    /*compute base level estimate*/
    dabcrs(abc_pl[0], data,thetal,rhol);
    for (j=0;j<M->G.numPoints;j++)
    {
        dmcint(abc_pl[0].nacc,abc_pl[0].k,thetal,M->G.coords + j*M->G.dim,
               M->g,E+j,V+j);
    }
    /* initial estimator*/
    for (j=0;j<M->G.numPoints;j++)
    {
        M->F[j] = E[j];
    }
    monotone(M,M_mon);
#if defined(__CHECKPOINT__)
        {
            FILE *fp;
            fp =fopen(CHECKPOINT_FILENAME,"a");
            /*save CDF state ... the Grid can be restored from main code*/
            fprintf(fp,"%d ",l);
            for (j=0;j<M->G.numPoints;j++)
            {
                fprintf(fp," %lg",M->F[j]);
            }
            fprintf(fp,"\n");
            /* save current empirical support*/
            fprintf(fp,"%d %d",l,abc_p.k);
            for (j=0;j<2*abc_p.k;j++)
            {
                fprintf(fp," %lg",abc_pl[l].support[j]);
            }
            fprintf(fp,"\n");
            /* save Nl sequence*/
            fprintf(fp,"%d %d",l, ml_p.L+1);
            for (j=0;j<=ml_p.L;j++)
            {
                fprintf(fp," %d",ml_p.Nl[j]);
            }
            fprintf(fp,"\n");
            /* save sim counter*/
            fprintf(fp,"%d %d\n",l,sim_counter);
            fclose(fp);
        }
#endif


    /*process bias updates*/
    for (l=1;l<=ml_p.L;l++)
    {
        memcpy(thetalm1,thetal,abc_pl[l-1].k*abc_pl[l-1].nacc*sizeof(double));
        memcpy(abc_pl[l].support,thetalm1,abc_pl[l].k*sizeof(double));
        memcpy(abc_pl[l].support+abc_pl[l].k,thetalm1,
               abc_pl[l].k*sizeof(double));
        /*update support for prior*/
        for (i=1;i<abc_pl[l-1].nacc;i++)
        {
            for (j=0;j<abc_pl[l-1].k;j++)
            {
                if (thetalm1[i*abc_pl[l].k +j] < abc_pl[l].support[j])
                {
                    abc_pl[l].support[j] = thetalm1[i*abc_pl[l].k + j];
                }

                if (thetalm1[i*abc_pl[l].k +j] > 
                    abc_pl[l].support[abc_pl[l].k + j])
                {
                    abc_pl[l].support[abc_pl[l].k + j] = 
                                                    thetalm1[i*abc_pl[l].k + j];
                }
            }
        }
        dabcrjs(abc_pl[l-1],abc_pl[l],data,M_mon, thetalm1,thetal,rhol);
        for (j=0;j<M->G.numPoints;j++)
        {
            dmcintd(abc_pl[l].nacc,abc_pl[l].k,thetal,thetalm1,
                    M->G.coords+j*M->G.dim,M->g,E+j,V+j);
        }
        
        for (j=0;j<M->G.numPoints;j++)
        {
            M->F[j] += E[j];
        }
        monotone(M,M_mon);
#if defined(__CHECKPOINT__)
        {
            FILE *fp;
            fp =fopen(CHECKPOINT_FILENAME,"a");
            /*save CDF state ... the Grid can be restored from main code*/
            fprintf(fp,"%d ",l);
            for (j=0;j<M->G.numPoints;j++)
            {
                fprintf(fp," %lg",M->F[j]);
            }
            fprintf(fp,"\n");
            /* save current empirical support*/
            fprintf(fp,"%d %d",l,abc_p.k);
            for (j=0;j<2*abc_p.k;j++)
            {
                fprintf(fp," %lg",abc_pl[l].support[j]);
            }
            fprintf(fp,"\n");
            /* save Nl sequence*/
            fprintf(fp,"%d %d",l, ml_p.L+1);
            for (j=0;j<=ml_p.L;j++)
            {
                fprintf(fp," %d",ml_p.Nl[j]);
            }
            fprintf(fp,"\n");
            /* save sim counter*/
            fprintf(fp,"%d %d\n",l,sim_counter);
            fclose(fp);
        }
#endif

    }
    
    memcpy(M->F,M_mon->F,M->G.numPoints*sizeof(double));

    free(thetal);
    free(thetalm1);
    free(rhol);
    free(abc_pl);
    free(E);
    free(V);
    return 0;
}



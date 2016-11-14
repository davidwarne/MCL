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
/**
 * @brief Approximate Bayesian computation partial rejection control.
 * @details Likelihood free sequential Monte Carlo as proposed by Sisson et al. [1]
 *
 * @param abc_p ABC Parameters
 * @param smc_p SMC parameters
 * @param data dataset to condition on
 * @param theta array to store resulting particles
 * @param weights normalised paticle weights
 * @param rho sample discrepency metric values
 *
 * @note [1] Sisson, S., Fan, Y. & Tanaka, M. M. 
 *          Sequential Monte Carlo without likelihoods.
 *          Proceedings of the National Academy of Sciences of the United States of America.
 *          2007;104(6):1760--1765.
 */
int dabcpcr(ABC_Parameters abc_p, SMC_Parameters smc_p, Dataset * data, double * theta, double * weights, double *rho)
{
    unsigned int t,i,j;
    Dataset *data_s;
    double *theta_prev;
    double * weights_prev;
    double * rho_prev;
    double W_sum, W_sum_prev;
    t =0;
    /*allocate memory for simulated  dataset*/
    data_s = copyDataset(data); 
    
    theta_prev = (double *)malloc(abc_p.k*abc_p.nacc*sizeof(double));
    weights_prev = (double *)malloc(abc_p.nacc*sizeof(double));
    rho_prev = (double *)malloc(abc_p.nacc*sizeof(double));

    /*initialise with ABC rejection samples with eps_0*/
    abc_p.eps = smc_p.eps_t[0];
    dabcrs(abc_p,data, theta, rho);
   
   /*initialise weights W_i = 1 (we don't normalise, rather we just store the sum)*/
    for (i=0;i<abc_p.nacc;i++)
    {
        weights[i] = 1.0;
    }
    W_sum = (double)abc_p.nacc;
#if defined(__CHECKPOINT__)
        {
            /*for long running simulations*/
            FILE *fp;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            /*output particles*/
            fprintf(fp,"%d %lg",t,smc_p.eps_t[t]);
            for (i=0;i<abc_p.nacc;i++)
            {
                for (j=0;j<abc_p.k;j++)
                {
                    fprintf(fp," %lg",theta[i*abc_p.k + j]);
                }
            }
            fprintf(fp,"\n");
            /*output weights*/
            fprintf(fp,"%d %lg",t,smc_p.eps_t[t]);
            for (i=0;i<abc_p.nacc;i++)
            {
                fprintf(fp," %lg",weights[i]);
            }
            fprintf(fp,"\n");

            fprintf(fp,"%d %d\n",t,sim_counter);
            fclose(fp);
        }
#endif
 
    /*commence sequential Monte Carlo steps*/
    for (t=1;t<smc_p.T;t++)
    {
        double ESS;
        memcpy(theta_prev,theta,abc_p.k*abc_p.nacc*sizeof(double));
        memcpy(weights_prev,weights,abc_p.nacc*sizeof(double));
        if (rho != NULL)
            memcpy(rho_prev,rho,abc_p.nacc*sizeof(double));
        W_sum_prev = W_sum;

        for (i=0;i<abc_p.nacc;i++)
        {
            double d, back_kern;
            d = INFINITY;
            while (d >= smc_p.eps_t[t])
            {
                /*sample a particle by weight from {theta_t-1,W_t-1} */
                j = durngpmfs(abc_p.nacc,weights_prev,W_sum_prev);
                /*perturb particle using transition kernel*/
                (*(smc_p.q))(abc_p.k,theta_prev+j*abc_p.k,theta + i*abc_p.k);
                /*simulate data*/
                (*(abc_p.s))(abc_p.sim,theta +i*abc_p.k,data_s);

                /*compute discrepency metric*/
                d = (*(abc_p.rho))(data,data_s);
            }

            if (rho != NULL)
                rho[i] = d;

            /*update particle weight*/
            back_kern = 0;
            for (j=0;j<abc_p.nacc;j++)
            {
                back_kern += weights_prev[j] *( (*(smc_p.qd))(abc_p.k,theta_prev+j*abc_p.k,theta+i*abc_p.k));
            }
            weights[i] = (*(abc_p.pd))(abc_p.k,theta + i*abc_p.k)*W_sum / back_kern;
            //{
            //    /*debug...*/
            //    FILE *fp;
            //    fp = fopen("")
            //}
            //fprintf(stderr,"%g,%g,%g,%g",weights_prev[i],weights[i],back_kern,(*(abc_p.pd))(abc_p.k,theta + i*abc_p.k));
            //for (j=0;j<abc_p.k;j++)
            //{
            //    fprintf(stderr,",%g,%g",theta_prev[i*abc_p.k + j],theta[i*abc_p.k + j]);
            //}
            //fprintf(stderr,"\n");
        }
        /*update weight sum*/
        W_sum = 0;
        for (i=0;i<abc_p.nacc;i++)
        {
            W_sum += weights[i];
        }
        
        /*compute effective sample size*/
        ESS = 0;
        for (i=0;i<abc_p.nacc;i++)
        {
            ESS += weights[i]*weights[i];
        }
        ESS = W_sum*W_sum/ESS;
        if (ESS < smc_p.E) /*resample with replacement*/
        {
            memcpy(theta_prev,theta,abc_p.k*abc_p.nacc*sizeof(double));
            if (rho != NULL)
                memcpy(rho_prev,rho,abc_p.nacc*sizeof(double));
            for (i=0;i<abc_p.nacc;i++)
            {
                j = durngpmfs(abc_p.nacc,weights,W_sum);
                memcpy(theta +i*abc_p.k, theta_prev + j*abc_p.k,abc_p.k*sizeof(double));
                if (rho != NULL)
                    rho[i] = rho[j];
            }
            /*reset weights*/
            for (i=0;i<abc_p.nacc;i++)
            {
                weights[i] = 1.0;
            }
            W_sum = (double)abc_p.nacc;
        }
#if defined(__CHECKPOINT__)
        {
            /*for long running simulations*/
            FILE *fp;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            /*output particles*/
            fprintf(fp,"%d %lg",t,smc_p.eps_t[t]);
            for (i=0;i<abc_p.nacc;i++)
            {
                for (j=0;j<abc_p.k;j++)
                {
                    fprintf(fp," %lg",theta[i*abc_p.k + j]);
                }
            }
            fprintf(fp,"\n");
            /*output weights*/
            fprintf(fp,"%d %lg",t,smc_p.eps_t[t]);
            for (i=0;i<abc_p.nacc;i++)
            {
                fprintf(fp," %lg",weights[i]);
            }
            fprintf(fp,"\n");

            fprintf(fp,"%d %d\n",t,sim_counter);
            fclose(fp);
        }
#endif
    }

    /*normalise weights*/
    for (i=0;i<abc_p.nacc;i++)
    {
        weights[i] /= W_sum;
    }
    return 0;
}

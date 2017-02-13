/* MCL: Monte Carlo Library 
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
 * @brief Pseudo-Marginal Metropolis-Hastings.
 * @details A standard Likelihood-free version of the Metropolis-Hastings
 * MCMC method. In the Pseudo-Marginal approach the true likelihood 
 * function p(X|theta) is replaced by an unbiased Monte Carlo estimator 
 * [1].
 *
 * @param pm_p Pseudo-Marginal parameters (i.e. provides estimator to 
 *             likelihood)
 * @param mcmc_p MCMC Parameters
 * @param data dataset to condition on
 * @param theta array to store Markov Chain state.
 *
 * @note [1] Andrieu, C. and Roberts, G. O. 
 *           The pseudo-marginal approach for efficient Monte Carlo 
 *           computations.
 *           The Annals of Statistics.
 *           2009;37(2):697--725.
 * @note For numerical stability reasons, all density computations are 
 * assumed to be in log scale.
 */
int
dmmh(PM_Parameters pm_p, MCMC_Parameters mcmc_p, Dataset *data, 
     double *theta)
{
    unsigned int cur_p, prop_p;
    double cur_logl, prop_logl;
    double cur_logp, prop_logp;
    double cur_logq, prop_logq;
    double log_h;
    /*initialise the chain*/
    memcpy(theta,mcmc_p.theta0,mcmc_p.k*sizeof(double));
    /*compute log-likelihood using Monte Carlo using theta'*/
    cur_logl = (*(pm_p.L_hat))(pm_p.alg_params,theta,data);
    cur_logp = (*(pm_p.pd))(mcmc_p.k,theta);
    /*perform burn-in iters*/
    cur_p = 0;
    prop_p = 1;
    for (int j=0;j<=mcmc_p.burnin_iters;j++)
    {
        /*sample proposal density q(theta -> theta')*/
        (*(mcmc_p.q))(mcmc_p.k,theta+cur_p*mcmc_p.k,theta+prop_p*mcmc_p.k);
        /*compute log-likelihood using Monte Carlo using theta'*/
        prop_logl = (*(pm_p.L_hat))(pm_p.alg_params,theta+prop_p*mcmc_p.k,
                                    data);
        /*Update log-prior*/
        prop_logp = (*(pm_p.pd))(mcmc_p.k,theta+prop_p*mcmc_p.k);
        /*update log-proposals*/
        prop_logq = (*(mcmc_p.qd))(mcmc_p.k,theta+prop_p*mcmc_p.k,
                                            theta+cur_p*mcmc_p.k);
        cur_logq = (*(mcmc_p.qd))(mcmc_p.k,theta+cur_p*mcmc_p.k,
                                           theta+prop_p*mcmc_p.k);
        /*compute log of MH acceptance probability*/
        log_h = prop_logl + prop_logp + prop_logq - cur_logl - cur_logp 
                - cur_logq;
        /*accept/reject*/
        if (log(durngus(0.0,1.0)) <= log_h )
        {
            unsigned int temp;
            /*since we are int burn-in we just swap indices*/
            temp = cur_p;
            cur_p = prop_p;
            prop_p = cur_p;

            /*update the cur_logl and cur_logp*/
            cur_logl = prop_logl;
            cur_logp = prop_logp;
        }
    }
    /*now start sampling*/
    cur_p = 0;
    prop_p = 1;
    for (int j=1;j<mcmc_p.iters;j++)
    {
        /*sample proposal density q(theta -> theta')*/
        (*(mcmc_p.q))(mcmc_p.k,theta+cur_p*mcmc_p.k,theta+prop_p*mcmc_p.k);
        /*compute log-likelihood using Monte Carlo using theta'*/
        prop_logl = (*(pm_p.L_hat))(pm_p.alg_params,theta+prop_p*mcmc_p.k,
                                    data);
        /*Update log-prior*/
        prop_logp = (*(pm_p.pd))(mcmc_p.k,theta+prop_p*mcmc_p.k);
        /*update log-proposals*/
        prop_logq = (*(mcmc_p.qd))(mcmc_p.k,theta+prop_p*mcmc_p.k,
                                            theta+cur_p*mcmc_p.k);
        cur_logq = (*(mcmc_p.qd))(mcmc_p.k,theta+cur_p*mcmc_p.k,
                                           theta+prop_p*mcmc_p.k);
        /*compute log of MH acceptance probability*/
        log_h = prop_logl + prop_logp + prop_logq - cur_logl - cur_logp 
                - cur_logq;

        /*accept/reject*/
        if (log(durngus(0.0,1.0)) <= log_h )
        {
            cur_logl = prop_logl;
            cur_logp = prop_logp;
        }
        else
        {
            /*copy previos state*/
            memcpy(theta+prop_p*mcmc_p.k,theta+cur_p*mcmc_p.k,
                   mcmc_p.k*sizeof(double));
        }
        cur_p = prop_p;
        prop_p++;
#if defined(__CHECKPOINT__)
        if(j%1000 == 0)
        {
            FILE *fp;
            unsigned int ii;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            for (int i=0;i<1000;i++)
            {
                ii = j - 1000 + i;
                fprintf(fp,"%d",ii);
                for (int k=0;k<abc_p.k;k++)
                {
                    fprintf(fp," %lg",theta[ii*mcmc_p.k+k]);
                }
                fprintf(fp,"\n");
            }
            fclose(fp);
        }
#endif
    }
    return 0;
}

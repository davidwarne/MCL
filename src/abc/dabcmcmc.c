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
 * @brief Approximate Bayesian computation Markov chain Monte Carlo.
 * @details implements the likelihood-free MCMC method of Marjoram et al. [1]
 *
 * @param abc_p ABC parameters
 * @param mcmc_p MCMC parameters
 * @param data dataset to condition on
 * @param theta array to store Markov Chain state (shoul represent samples from the posterior)
 * @param rho sample discrepency metric values.
 *
 * @note [1] Marjoram, P., Molitor, J., Plagnol, V. & Tavare, S. 
 *          Markov chain Monte Carlo without likelihoods.  
 *          Proceedings of the National Academy of Sciences of the United States of America. 
 *          2003;100(26):15324--15328.   
 */
int dabcmcmc(ABC_Parameters abc_p, MCMC_Parameters mcmc_p, Dataset *data, double * theta, double *rho)
{
    unsigned int j;
    unsigned int cur_p,prop_p;
    unsigned int nacc;    
    double d;
    Dataset *data_s;
    /*allocate memory for simulated dataset*/
    data_s = copyDataset(data);
    
    
    if (mcmc_p.theta0 == NULL)
    {
        /*obtain initial accepted sample using ABC rejection*/
        nacc = abc_p.nacc;
        abc_p.nacc = 1;
        dabcrs(abc_p,data,theta,rho);
        abc_p.nacc = nacc;
    }
    else
    {
        memcpy(theta,mcmc_p.theta0,abc_p.k*sizeof(double));
    }


    /*now simulate the markov chain*/
    cur_p = 0;
    prop_p = 1;
    for (j = 0; j<=mcmc_p.burnin_iters;j++)
    {
        /*make new proposal*/
        (*(mcmc_p.q))(abc_p.k,theta+cur_p*abc_p.k,theta+prop_p*abc_p.k);

        /*simulate process with proposed parameters*/
        (*(abc_p.s))(abc_p.sim,theta+prop_p*abc_p.k,data_s);

        /*check if we accept the new proposal*/
        d = (*(abc_p.rho))(data,data_s);
        if (d <= abc_p.eps)
        {
            double h,lr,p_n,q_n,p_d,q_d;
            
            /*prior prob of proposal*/ 
            p_n = (*(abc_p.pd))(abc_p.k,theta + prop_p*abc_p.k);             
            
            /*prior prob of current state*/ 
            p_d = (*(abc_p.pd))(abc_p.k,theta + cur_p*abc_p.k);             
            
            /*backward transition probability*/
            q_n = (*(mcmc_p.qd))(abc_p.k,theta+prop_p*abc_p.k,theta+cur_p*abc_p.k);
            
            /*forward transition probability*/
            q_d = (*(mcmc_p.qd))(abc_p.k,theta+cur_p*abc_p.k,theta+prop_p*abc_p.k);
            
            /*p(th')q(th' -> th)/p(th)q(th -> th') likelihood ratio*/
            lr = (p_n*q_n)/(p_d*q_d);
            h = (lr < 1.0) ? lr : 1.0;

            if (durngus(0.0,1.0) <= h) /*accept proposal*/
            {
                unsigned int temp;
                temp = cur_p; /*we are in burn-in phase, so just swap indices*/
                cur_p = prop_p;
                prop_p = cur_p;
            }
        }
    }
    cur_p = 0;
    prop_p = 1;
    for (j=1;j<abc_p.nacc;j++)
    {
        unsigned char reject;
        reject = 1;
        /*make new proposal*/
        (*(mcmc_p.q))(abc_p.k,theta+cur_p*abc_p.k,theta+prop_p*abc_p.k);

        /*simulate process with proposed parameters*/
        (*(abc_p.s))(abc_p.sim,theta+prop_p*abc_p.k,data_s);

        /*check if we accept the new proposal*/
        d = (*(abc_p.rho))(data,data_s);
        if (d <= abc_p.eps)
        {
            double h,lr,p_n,q_n,p_d,q_d;
            
            /*prior prob of proposal*/ 
            p_n = (*(abc_p.pd))(abc_p.k,theta + prop_p*abc_p.k);             
            
            /*prior prob of current state*/ 
            p_d = (*(abc_p.pd))(abc_p.k,theta + cur_p*abc_p.k);             
            
            /*backward transition probability*/
            q_n = (*(mcmc_p.qd))(abc_p.k,theta+prop_p*abc_p.k,theta+cur_p*abc_p.k);
            
            /*forward transition probability*/
            q_d = (*(mcmc_p.qd))(abc_p.k,theta+cur_p*abc_p.k,theta+prop_p*abc_p.k);
            
            /*p(th')q(th' -> th)/p(th)q(th -> th') likelihood ratio*/
            lr = (p_n*q_n)/(p_d*q_d);
            h = (lr < 1.0) ? lr : 1.0;

            reject = (durngus(0.0,1.0) > h); /*reject condition*/
        }
        
        if (reject)/*not a candidate, keep current state*/
        {
            memcpy(theta+prop_p*abc_p.k,theta+cur_p*abc_p.k,abc_p.k*sizeof(double));
            if (rho != NULL)
                rho[prop_p] = rho[cur_p];
        }
        else /*accept new state*/
        {
            if (rho != NULL)
                rho[prop_p] = d;
        }
        /*update pointers*/
        cur_p = prop_p;
        prop_p++;
#if defined(__CHECKPOINT__)
        if(j%1000 == 0)
        {
            
            FILE *fp;
            unsigned int i,ii,k;
            fp = fopen(CHECKPOINT_FILENAME,"a");
            for (i=0;i<1000;i++)
            {
                ii = j - 1000 + i;
                fprintf(fp,"%d",ii);
                for (k=0;k<abc_p.k;k++)
                {
                    fprintf(fp," %lg",theta[ii*abc_p.k+k]);
                }
                fprintf(fp,"\n");
            }
            fclose(fp);
        }
#endif
    }

    return 0;
}

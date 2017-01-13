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
 * @brief ABC rejection sampling of joint distribution
 * @deials approximation to strong converging sequence as per Warne et al. 
 * used of MLABC bias correction step. 
 *
 * @param abc_p_c ABC parameters for coarse level
 * @param abc_p_f ABC parameters for fine level
 * @param data dataset X_d
 * @param theta_c array of accepted posterior samples coupled with theta_f
 * @param theta_f array of accepted posterior samples at fine level
 * @param rho_f values of rho
 */
int 
dabcrjs(ABC_Parameters abc_p_c, ABC_Parameters abc_p_f,Dataset *data, 
        CDF_estimate *F_c, double *theta_c, double *theta_f, double * rho_f)
{
    size_t *G_LUT; /*grid lookup*/
    size_t ind_f;
    size_t ind_c;
    size_t sub_f[255]; 
    size_t sub_c[255]; 
    double *F_tmp;
    unsigned int i,j,k;
    /* use standard ABC rejection for fine level sample*/
    dabcrs(abc_p_f,data,theta_f,rho_f);
    
    G_LUT = (size_t *)malloc(abc_p_f.nacc*sizeof(size_t));

    F_tmp = (double*)malloc(F_c->G.numPoints*sizeof(double));

    /*compute empirical CDF and identify samples with grid cell*/
    for (j=0;j<F_c->G.numPoints;j++)
    {
        F_tmp[j] = 0;
        for (i=0;i<abc_p_f.nacc;i++)
        {
            unsigned int tf;
            double s[abc_p_f.k];
            tf = 1;
            for (k=0;k<abc_p_f.k;k++)
            {
                tf = tf && 
                   (theta_f[i*abc_p_f.k + k] <= F_c->G.coords[j*abc_p_f.k + k]);
                tf = tf && 
                    (theta_f[i*abc_p_f.k + k] > 
                           (F_c->G.coords[j*abc_p_f.k + k] - F_c->G.deltas[k]));
            }
            F_tmp[j] += (*(F_c->g))(abc_p_f.k,theta_f +i*abc_p_f.k,
                                    F_c->G.coords + j*F_c->G.dim);
            if (tf == 1)
            {
                G_LUT[i] = j;
            }
        }
    }
    for (j=0;j<F_c->G.numPoints;j++)
    {
        F_tmp[j] /= ((double)abc_p_f.nacc);
    }

    /*match marginals to determine coarse samples*/
    for (i=0;i<abc_p_f.nacc;i++)
    {
        /*convert index to grid coordinate*/
        ind_f = G_LUT[i];
        for (k=0;k<abc_p_f.k;k++)
        {
            sub_f[k] = ind_f % F_c->G.D;
            ind_f = ind_f / F_c->G.D;
        }
        /*locate grid point in coarse level estimator closest wrt marginal 
         * probabilities*/
        for (k=0;k<abc_p_f.k;k++)
        {
            int d;
            double dist_min = INFINITY;
            for (d=0;d<F_c->G.D;d++)
            {
                double dist;
                
                dist = fabs(F_tmp[F_c->G.offsets[k] + sub_f[k]*F_c->G.incs[k]] 
                            - F_c->F[F_c->G.offsets[k] + d*F_c->G.incs[k]]);
                if (dist < dist_min)
                {
                    dist_min = dist;
                    sub_c[k] = d;
                }
            }
        }
        /*convert grid coordinate to index*/
        ind_c = 0;
        for (k=0;k<abc_p_c.k;k++)
        {
            ind_c += sub_c[k]*F_c->G.incs[k];
        }

        /*select this point + perturbation*/
        for( k=0;k<abc_p_c.k;k++)
        {
           theta_c[i*abc_p_c.k + k] = F_c->G.coords[ind_c*abc_p_c.k + k] 
                                         + durngus(-1.0,1.0)*(F_c->G.deltas[k]);
        }
    }
    return 0;
}


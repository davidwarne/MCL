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
 * @brief standard ABC rejection
 * @detail Applies ABC rejection to compute the expectation of an indicator 
 * function at various points in parameter space resulting in an estimator of 
 * the posterior CDF.
 *
 * @param sim SSAL simulation
 * @param abc_p ABC parameters
 * @param data a realisation of the system we a fitting to
 * @param M the resulting approximated CDF
 */
int 
dabccdfs(ABC_Parameters abc_p,  Dataset *data, CDF_estimate * M)
{
    double *theta;
    unsigned int j;

    /*samples allocation*/
    if ((theta = (double *)malloc(abc_p.nacc*abc_p.k*sizeof(double))) == NULL) 
    {
        return 1;
    }
    /*ABC rejection sampling*/
    dabcrs(abc_p, data,theta, NULL);
    /*collect into indicator functions*/
    for (j=0;j<M->G.numPoints;j++)
    {
        dmcint(abc_p.nacc,abc_p.k,theta,M->G.coords +j*M->G.dim,M->g,
               &(M->F[j]),&(M->V[j]));
    }
    return 0;
}


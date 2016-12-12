/* ABC: approximate Bayesian Computation
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

#ifndef __LIBABC_H_
#define __LIBABC_H_

#include "SSAL.h"
#include <math.h>
#include <stdint.h>

#define ABC_MAX_NAME_SIZE 128
enum ABC_DataType_enum{
    REAL32_DATA,
    REAL64_DATA,
    INT32_DATA,
    INT64_DATA,
    NAT32_DATA,
    NAT64_DATA,
    TEXT_DATA
};
typedef enum ABC_DataType_enum ABC_DataType;

#if defined(__CHECKPOINT__)
#define CHECKPOINT_FILENAME "sim.chkpnt"
#endif
extern unsigned int sim_counter; 

struct field_struct {
    char name[ABC_MAX_NAME_SIZE]; /*field name*/
    size_t numRows; 
    size_t numCols;
    size_t numBytes;
    ABC_DataType type;
    void *data_array;

};
typedef struct field_struct field;

/* generic dataset structure*/
struct Dataset_struct {
    size_t numFields;
    field *fields;
};
typedef struct Dataset_struct Dataset;

struct ABC_Parameters_struct {
    unsigned int nmax; /*total max samples before termination*/
    SSAL_real_t percentile; /*If nmax is used, then use percentile to INFER eps*/
    SSAL_real_t eps; /*acceptance threshold*/
    unsigned int nacc; /*number of acceptances required*/
    unsigned int k; /*dimensionality of parameter space*/
    SSAL_real_t *support; /*prior sampling region*/
    void *sim; /*Model simulation object */
    SSAL_real_t (*rho)(Dataset *,Dataset *); /*distance function*/
    int (*p)(unsigned int,unsigned int,SSAL_real_t *,SSAL_real_t *); /*prior sampler function*/
    SSAL_real_t (*pd)(unsigned int, SSAL_real_t*); /* Prior Distribution function; not required for all methods*/
    int (*s)(void *, SSAL_real_t* ,Dataset *); /*model simulation function*/
};

typedef struct ABC_Parameters_struct ABC_Parameters;

struct MCMC_Parameters_struct {
    unsigned int burnin_iters; /*number of burnin simulation to perfrom before sampling*/
    SSAL_real_t *theta0;
    int (*q)(unsigned int, SSAL_real_t*, SSAL_real_t*); /*transition kernel sampler*/
    SSAL_real_t (*qd)(unsigned int, SSAL_real_t*, SSAL_real_t*); /*transition kernel Distribution function*/
};
typedef struct MCMC_Parameters_struct MCMC_Parameters;

struct SMC_Parameters_struct {
    unsigned int T;/*number of steps*/
    SSAL_real_t *eps_t; /*ABC thresholds*/
    SSAL_real_t E; /*Effective sample size (ESS) threshold for re-sampling.*/
    int (*q)(unsigned int, SSAL_real_t*, SSAL_real_t*); /*transition kernel sampler*/
    SSAL_real_t (*qd)(unsigned int, SSAL_real_t*, SSAL_real_t*); /*transition kernel Distribution function*/
};
typedef struct SMC_Parameters_struct SMC_Parameters;

struct MLMC_Parameters_struct {
    unsigned int L; /*Max level*/
    unsigned int K; /* scale factor */
    SSAL_real_t eps0; /*coarsest level*/
    SSAL_real_t *eps_l; /*sequence of thresholds*/
    unsigned int *Nl; /*sequence of Monte Carlo sample numbers per level*/
    SSAL_real_t target_RMSE; /*target error (in root mean square)*/
    SSAL_real_t gamma; /*coupled sampler complexity order*/
    SSAL_real_t alpha; /*weak convergence rate*/
    SSAL_real_t beta; /*strong convergence rate*/
    unsigned int presample; /* Flag to use sampling for variance estimates*/
    unsigned int presample_trials; /*number of samples for variance estimates*/
};
typedef struct MLMC_Parameters_struct MLMC_Parameters;

struct ndgrid_struct {
    int dim; /*dimensionality of grid*/
    int D; /*points per axis*/
    unsigned int numPoints; /*total grid points D^dim*/
    double *coords;
    double deltas[255];
    /*axis offset abd incerements*/
    size_t offsets[255];
    size_t incs[255];
};
typedef struct ndgrid_struct ndgrid;

struct CDF_estimate_struct {
    ndgrid G; /*nd-grid*/
    double *F; /*probabilities F[j] = E[g((theta -x)/delta)]*/
    double *V; /*variance V[g((theta -x)/delta)]*/
    double (*g)(int,double*,double*); /*smoothing functional*/
};
typedef struct CDF_estimate_struct CDF_estimate;

#define SUB2IND(sub,dim,D,ind) {\
    int k;\
    (ind) = 0;\
    for (k=(dim)-1;k>=0;k--) {\
        (ind) *= (D);\
        (ind) += (sub)[k];\
    }\
}

#define IND2SUB(ind,dim,D,sub) {\
    int k, ind_tmp; \
    ind_tmp = (ind); \
    for (k=0;k<(dim);k++) {\
        (sub)[k] = ind_tmp % (D); \
        ind_tmp = ind_tmp / (D); \
    }\
}


/*ABC CDF estimators*/
//int sabccdfs(SSAL_Simulation *,ABC_Parameters,  float *, int, float *,float *,float *);

//int smlabccdfs(SSAL_Simulation *,ABC_Parameters ,ML_Parameters,  float *, int, float *,float *,float *,unsigned int *, float **);

int dabccdfs(ABC_Parameters, Dataset * , CDF_estimate *);

int dmlabccdfs(ABC_Parameters ,MLMC_Parameters,  Dataset *, CDF_estimate *);

/*ABC rejection samplers*/
int dabcrs(ABC_Parameters, Dataset *, double *,double *);
int dabcrjs(ABC_Parameters , ABC_Parameters, Dataset * , CDF_estimate *, double *, double *, double *);

int dabcmcmc(ABC_Parameters, MCMC_Parameters, Dataset *, double *,double *);
int dabcpcr(ABC_Parameters, SMC_Parameters, Dataset *, double *, double *,double *);

/*sampling calculation of N_l*/
int dmlabcnls(ABC_Parameters ,int ,MLMC_Parameters,  Dataset *, CDF_estimate *,unsigned int *,double *);
int dabcns(ABC_Parameters, int, double ,Dataset * , CDF_estimate *, unsigned int *);

/*Monte Carlo integration functions*/
int dmcint(unsigned int, unsigned int, double *, double*, double (*)(int , double *,double *), double *, double *); 
int dmcintd(unsigned int, unsigned int, double *, double *, double *,double (*)(int , double *, double*), double *, double *); 

Dataset *copyDataset(Dataset *);
void monotone(CDF_estimate *, CDF_estimate *);
double cdf2pmf(CDF_estimate *, double *);
int dmlabcnlps(ABC_Parameters abc_p, int trials,MLMC_Parameters ml_p, Dataset *data, CDF_estimate *M, double *Yl, double *Ylu, double *Fl);
#endif

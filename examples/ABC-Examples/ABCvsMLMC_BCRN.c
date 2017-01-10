/* LIBABC: approximate Bayesian Computation
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

#include "mcl.h"
#include <time.h>


/*prior bounds*/
SSAL_real_t *a;
SSAL_real_t *b;

/*smoothing deltas*/
SSAL_real_t *deltas;
/*data set*/
Dataset dataset;

/*filenames*/
char *model_filename;
char *data_filename;
char *outputfilename;
char *summaryfilename;
/*parameter structures*/
ABC_Parameters abc_p;
MLMC_Parameters ml_p;
SSAL_real_t smoothfactor;
SSAL_real_t scalefactor;
int D; 
char ml;

/* function prototypes*/
Dataset ImportData(char *,int);
SSAL_real_t rho(Dataset *, Dataset *);
int prior(unsigned int, unsigned int,SSAL_real_t *, SSAL_real_t *);
int simulate_model(void *, SSAL_real_t *, Dataset *);
ndgrid GenNdGrid(int, SSAL_real_t *,SSAL_real_t *,int,double *);
SSAL_real_t g(int , SSAL_real_t *,SSAL_real_t *);
SSAL_real_t I(int , SSAL_real_t *,SSAL_real_t *);
int getCovergenceParams(char *,unsigned int, SSAL_real_t *, SSAL_real_t *);
void GetCmdLineArgs(int,char **);

/*program entry point*/
int main(int argc, char ** argv)
{
    SSAL_Model model;
    SSAL_ChemicalReactionNetwork *CRN_ptr;
    SSAL_RealisationSimulation *RS_ptr;
    SSAL_Simulation sim;
    CDF_estimate M;
    double *T;
    double *X_d;
    int i,j,l;
    clock_t start_t,end_t;
    clock_t M_t,ML_t;
    
    /*some defaults*/
    smoothfactor = 1.0;
    scalefactor = 10.0;
    ml = 0;
    ml_p.presample = 1;
    
    /*get user args*/
    GetCmdLineArgs(argc,argv); 
    /*init SSAL */
    SSAL_Initialise(argc,argv);
    
    /*read model defintion*/
    model = SSAL_ImportLSBML(model_filename);
    CRN_ptr = (SSAL_ChemicalReactionNetwork *)model.model;
    SSAL_WriteChemicalReactionNetwork(stdout,*CRN_ptr);
    
    /*read data*/
    dataset = ImportData(data_filename,CRN_ptr->N);
    T = (double*)dataset.fields[0].data_array;
    X_d = (double*)dataset.fields[1].data_array;

    /*allocate prior dist bounds*/
    a = (SSAL_real_t *)malloc(CRN_ptr->M*sizeof(SSAL_real_t));
    b = (SSAL_real_t *)malloc(CRN_ptr->M*sizeof(SSAL_real_t));
    deltas = (SSAL_real_t *)malloc(CRN_ptr->M*sizeof(SSAL_real_t));
    for (j=0;j<CRN_ptr->M;j++)
    {
        a[j] = 0;
        b[j] = (CRN_ptr->c[j])*scalefactor;
        deltas[j] = smoothfactor*(b[j] - a[j])/((SSAL_real_t)(D -1));
    }
    /*gamma = number of timesteps, J = number of indicator functions*/
    ml_p.gamma = dataset.fields[0].numRows*CRN_ptr->N;
    ml_p.eps_l = (SSAL_real_t*)malloc((ml_p.L+1)*sizeof(SSAL_real_t));
    /*init epsilon sequence*/
    ml_p.eps_l[0] =  ml_p.eps0;
    for (l=1;l<=ml_p.L;l++)
    {
        ml_p.eps_l[l] =  ml_p.eps0/powf((SSAL_real_t)(ml_p.K),(SSAL_real_t)l);
        //ml_p.eps_l[l] =  ml_p.eps0/powf((SSAL_real_t)(ml_p.K),(SSAL_real_t)l);
    }
    /*initialise equivalent ABC parameters*/
    abc_p.nmax = 0;
    abc_p.eps = ml_p.eps_l[ml_p.L]; 
    abc_p.k = CRN_ptr->M;
    abc_p.rho = &rho;
    abc_p.p = &prior;
    abc_p.s = &simulate_model;
    abc_p.sim = (void*)&sim;
    
    /*build simulation*/
    sim = SSAL_CreateRealisationsSim(&model,CRN_ptr->N,NULL,1,dataset.fields[0].numRows,T,CRN_ptr->X0);
    RS_ptr = (SSAL_RealisationSimulation *)(sim.sim);
    /*generate simulation data*/
    SSAL_Simulate(&sim,SSAL_ESSA_GILLESPIE_SEQUENTIAL,NULL);
    for (i=0;i<dataset.fields[0].numRows;i++)
    {
        fprintf(stdout,"X_d(%f) = ",T[i]);
        for (j=0;j<dataset.fields[1].numRows;j++)
        {
            fprintf(stdout," %f ",X_d[j*dataset.fields[1].numCols + i]);
            fprintf(stdout," %f ",RS_ptr->output[j*dataset.fields[1].numCols + i]);
        }
        fprintf(stdout,"\n");
    }

    /* Generate Grid and memory for CDF*/
    M.G = GenNdGrid(CRN_ptr->M,a,b,D,deltas);
    M.F = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.V = (double*)malloc(M.G.numPoints*sizeof(double)); 
    M.g = &g; 
    /*for support constraint update*/
    abc_p.support = (SSAL_real_t *)malloc(abc_p.k*2*sizeof(double));
    memcpy(abc_p.support,a,abc_p.k*sizeof(SSAL_real_t));
    memcpy(abc_p.support+abc_p.k,b,abc_p.k*sizeof(SSAL_real_t));

    if (ml)
    {
        int rc;
        FILE *fp;
        ml_p.Nl = (unsigned int *)malloc((ml_p.L+1)*sizeof(unsigned int));

        for (i=0;i<M.G.dim;i++)
        {
            fprintf(stdout,"d_%d = %g\n",i,M.G.deltas[i]);
        }
        /*trial samples to compute N_l*/
        dmlabcnls(abc_p,1000,ml_p,&dataset,&M,ml_p.Nl,NULL);
        
        for (l=0;l<=ml_p.L;l++)
        {
            fprintf(stdout,"N_l = %d,eps_l = %f\n",ml_p.Nl[l],ml_p.eps_l[l]);
        }
        /*run an time MLABC*/
        start_t = clock();
        rc = dmlabccdfs(abc_p,ml_p,&dataset,&M);
        end_t = clock();
        ML_t = end_t - start_t;
        
        /*ouput estimated joint CDF*/
        fp = fopen(outputfilename,"w");
        for (j=0;j<M.G.numPoints;j++)
        {
            fprintf(fp,"%d," ,j);
            for (i=0;i<M.G.dim;i++)
            {
                fprintf(fp,"%g,",M.G.coords[j*M.G.dim+i]);
            }
            fprintf(fp,"%g\n",M.F[j]);
        }
        fclose(fp);
       
        /*output summary data*/
        fp = fopen(summaryfilename,"w");
        fprintf(fp,"%g,%d,%u,%u,%g,%g,%g,%g,%g",((double)ML_t)/((double)CLOCKS_PER_SEC),M.G.numPoints,ml_p.K,ml_p.L,ml_p.eps0,ml_p.target_RMSE,ml_p.gamma,ml_p.alpha,ml_p.beta);
        fprintf(fp,"\n");
        fclose(fp);
    }
    else
    {
        int rc;
        FILE *fp;
        unsigned int N;
        
        if (ml_p.presample)
        {
            dabcns(abc_p,MLMC_TIMING_TRIALS,ml_p.target_RMSE,&dataset,&M,&N);
        }
        else
        {
            abc_p.nacc = ceil(log(M.G.numPoints)/(ml_p.target_RMSE*ml_p.target_RMSE));
        }
        
        fprintf(stdout,"N_ABC = %d\n",N);
        fprintf(stdout,"eps = %f\n",abc_p.eps);
        /*run and time standard ABC rejection*/
        start_t = clock();
        dabccdfs(abc_p,&dataset,&M);
        end_t = clock();
        M_t = end_t - start_t;
 
        /*ouput estimated joint CDF*/
        fp = fopen(outputfilename,"w");
        for (j=0;j<M.G.numPoints;j++)
        {
            fprintf(fp,"%d," ,j);
            for (i=0;i<M.G.dim;i++)
            {
                fprintf(fp,"%g,",M.G.coords[j*M.G.dim+i]);
            }
            fprintf(fp,"%g\n",M.F[j]);
        }
        fclose(fp);
       
        /*output summary data*/
        fp = fopen(summaryfilename,"w");
        fprintf(fp,"%g,%d,%g,%g\n",((double)M_t)/((double)CLOCKS_PER_SEC),M.G.numPoints,abc_p.eps,ml_p.target_RMSE);
        fclose(fp);
    }
}

/*import CRN realisation data*/
Dataset ImportData(char *filename,int n)
{
    FILE *fp;
    Dataset data;
    char c;
    int numlines;
    double* T;
    double *X;
    int t,i;
    fp = fopen(filename,"r");

    /*count lines*/
    numlines = 0;
    while (!feof(fp)) 
    {
        c = fgetc(fp);
        numlines += (c == '\n');
    }
    rewind(fp);

    data.numFields = 2; /*Field 1 is T, Field 2 is X(0:T)*/
    data.fields = (field *)malloc(data.numFields*sizeof(field));
    
    /*field 1: timesteps*/
    strncpy(data.fields[0].name,"T",ABC_MAX_NAME_SIZE);
    data.fields[0].numRows = numlines-2;
    data.fields[0].numCols = 1;
    data.fields[0].type = REAL64_DATA;
    data.fields[0].numBytes = 8*(numlines-2);
    data.fields[0].data_array = (void *)malloc(data.fields[0].numRows*data.fields[0].numCols*sizeof(double));
    
    /*field 1: timesteps*/
    strncpy(data.fields[1].name,"X_d",ABC_MAX_NAME_SIZE);
    data.fields[1].numRows = n;
    data.fields[1].numCols = numlines-2;
    data.fields[1].type = REAL64_DATA;
    data.fields[1].numBytes = 8*n*(numlines-2);
    data.fields[1].data_array = (void *)malloc(data.fields[1].numRows*data.fields[1].numCols*sizeof(double));

    /*read first line*/
    while ((c = fgetc(fp)) != '\n' ) ;
    while ((c = fgetc(fp)) != '\n' ) ;

    T = (double *)data.fields[0].data_array;
    X = (double *)data.fields[1].data_array;
    for (t=0;t<data.fields[0].numRows;t++)
    {
        fscanf(fp,"%lf",&(T[t]));
        for (i=0;i<data.fields[1].numRows;i++){
            fscanf(fp,",%lf",&(X[i*data.fields[1].numCols + t]));
        } 
        while ((c = fgetc(fp)) != '\n' ) ;
    }
    return data;
}

/*distance function*/
SSAL_real_t rho(Dataset *D, Dataset *D_s)
{
    int i,t;
    size_t n,nt;
    SSAL_real_t X_d_norm, X_diff_norm;
    SSAL_real_t *X_d, *X;
    SSAL_real_t d;
  
    n = D->fields[1].numCols;
    nt = D->fields[1].numRows;
    X_d = (SSAL_real_t *)D->fields[1].data_array;
    X = (SSAL_real_t *)D_s->fields[1].data_array;
    d = 0;
    for (t=0;t<nt;t++)
    {
        X_d_norm = 0;
        X_diff_norm = 0;
    //    fprintf(stdout,"\n ");
        for (i=0;i<n;i++)
        {
      //      fprintf(stdout,"%f %f\n",X_d[i*nt+t],X[i*nt+t]);
           X_d_norm += X_d[i*nt+t]*X_d[i*nt+t];
        }
        for (i=0;i<n;i++)
        {
            X_diff_norm += (X_d[i*nt + t] - X[i*nt+t])*(X_d[i*nt+t] - X[i*nt+t]);
        }

        if (X_d_norm != 0)
        {
            d += (X_diff_norm)/X_d_norm;
        }
        else
        {
            d += X_diff_norm;
        }
    }
    d /= (SSAL_real_t)nt;
       // fprintf(stdout,"%f\n ",sqrt(d));
    return sqrt(d);    
}

/*prior dist sampler*/
//int prior(unsigned int k, unsigned int numsamples, SSAL_real_t * theta)
int prior(unsigned int k, unsigned int numsamples,SSAL_real_t * support, SSAL_real_t * theta)
{
    int i,j;
    if (support == NULL)
  //  if (1)
    {
        for (i=0;i<numsamples;i++)
        {
            for (j=0;j<k;j++)
            {
                theta[k*i + j] = durngus(a[j],b[j]);
            }
        }
    } 
    else 
    {
        for (i=0;i<numsamples;i++)
        {
            for (j=0;j<k;j++)
            {
                theta[k*i + j] = durngus(support[j],support[k+j]);
            }
        }
    }
    return 0;
}

/*model simulation routine*/
int simulate_model(void *sim, SSAL_real_t * theta, Dataset * D_s)
{
    SSAL_real_t * X_s;
    SSAL_Simulation *sim_ptr;
    SSAL_RealisationSimulation *RS_ptr;
    SSAL_ChemicalReactionNetwork *CRN_ptr;

    /*cast pointers*/
    X_s = (SSAL_real_t *)D_s->fields[1].data_array;
    sim_ptr = (SSAL_Simulation *)sim;
    RS_ptr = (SSAL_RealisationSimulation *)(sim_ptr->sim);
    CRN_ptr = (SSAL_ChemicalReactionNetwork *)(sim_ptr->model->model);
    /*copy parameter vector to reaction rate vector*/
    memcpy(CRN_ptr->c,theta,CRN_ptr->M*sizeof(SSAL_real_t));
    /*Simulate with Gillespie direct method*/
    SSAL_Simulate(sim_ptr,SSAL_ESSA_GILLESPIE_SEQUENTIAL,NULL);
    memcpy(X_s,RS_ptr->output,RS_ptr->Nvar*RS_ptr->NT*sizeof(SSAL_real_t)); 
}


/*generate a regular N-D grid*/
ndgrid GenNdGrid(int d, SSAL_real_t *S_l,SSAL_real_t *S_u,int D,double * deltas)
{
    SSAL_real_t * axes;
    size_t *index;
    int i,j;
    ndgrid G; 

    index = (size_t *)malloc(d*sizeof(size_t));
    axes = (SSAL_real_t *)malloc(d*D*sizeof(SSAL_real_t));
    G.dim = d;
    G.D = D;
   
    memcpy(G.deltas,deltas,G.dim*sizeof(double));

    for (j=0;j<d;j++)
    {
        for (i=0;i<D;i++)
        {
            axes[j*D+i] = (((SSAL_real_t)i)*(S_u[j] - S_l[j]))/((SSAL_real_t)(D-1)) + S_l[j];
        }
    }

    G.numPoints = 1;
    for (j=0;j<d;j++)
    {
        G.numPoints *= D;
    }

    G.coords = (SSAL_real_t *)malloc(G.numPoints*d*sizeof(SSAL_real_t));
    /*init index*/
    for (j=0;j<d;j++)
    {
        index[j] = 0;
        G.coords[j] = axes[j*D + index[j]];
    }
    
    int c; /*carry flag*/
    /*generate grid points*/
    for (i=1;i<G.numPoints;i++)
    {
        c = 1;
        for (j=0;j<d;j++) /*base D adder with carry*/
        {
            index[j] = (index[j] + c) % D;
            c = (index[j] == 0 && c != 0);
        }
        /* generate point*/

        for (j=0;j<d;j++)
        {
            G.coords[i*d + j] = axes[j*D + index[j]];
        }
    }

    /*compute offsets and increments (for marginals) */
    for (j=0;j<d;j++)
    {
        int mult;
        mult = 1;
        G.offsets[j] = 0;
        for (i=0;i<d;i++)
        {
            if (i == j)
            {
                G.incs[j] = mult;
            }
            else
            {
                G.offsets[j] += (D-1)*mult;
            }
            mult *= D;
        }
    }

    free(index);
    free(axes);
    return G;
    
}

/*smoothing function*/
SSAL_real_t g(int d, SSAL_real_t *y, SSAL_real_t *x)
{
    int j;
    SSAL_real_t res = 1;
    SSAL_real_t s[255];

    for (j=0;j<d;j++)
    {
        s[j] = (y[j] -x[j])/deltas[j];    
        //s[j] = (y[j] -x[j]);    
        /*based on Giles method*/  
        res *= (s[j] <= 1) ? ((s[j] >= -1) ? 0.625*s[j]*s[j]*s[j] - 1.125*s[j] + 0.5 : 1 ) : 0;
    }
    return res;
}

SSAL_real_t gbs(int d, SSAL_real_t *y,SSAL_real_t * x)
{
    int j;
    SSAL_real_t res = 1;
    SSAL_real_t s[255];

    for (j=0;j<d;j++)
    {
        s[j] = (y[j] -x[j])/deltas[j];    
        //s[j] = (y[j] -x[j]);    
        /*cubic spline g' = 0 at edges*/  
        res *= (s[j] <= 1) ? ((s[j] >= -1) ? (s[j]+1.0)*(s[j]+1.0)*(s[j] + 1.0)/4.0 - (3.0*(s[j] + 1.0)*(s[j] + 1.0))/4.0 + 1.0 : 1 ) : 0;
    }
    return res;
}


/*smoothing function*/
SSAL_real_t I(int d, SSAL_real_t *y, SSAL_real_t *x)
{
    int j;
    SSAL_real_t res = 1;
    SSAL_real_t s[255];

    unsigned char tf;
    tf = 1;
    for (j=0;j<d;j++)
    {
      //  fprintf(stdout,"%f",x[j]);
        tf = tf && (y[j] <= x[j]);    
    }
    //fprintf(stdout,"\n");
    return (double)tf;
}

SSAL_real_t l_inf(int d, SSAL_real_t *y, SSAL_real_t *x)
{
    int j;
    SSAL_real_t n;
    SSAL_real_t en;
    n = 0;
    for (j=0;j<d;j++)
    {
        en = fabs(y[j]);
        n = (n <  en ) ? en : n; 
    }
    return n;
}


/**
 * Reads parameter file, and gets the empirical parameters for the data dimensionality
 */
int getCovergenceParams(char * filename,unsigned int d, SSAL_real_t *alpha, SSAL_real_t *beta)
{
    FILE *fp;
    unsigned int i;
    float alphaf,betaf;
    unsigned int df;
    df = 0;
    fp = fopen(filename,"r");

    while (d != df)
    {
        fscanf(fp,"%u,%f,%f\n",&df,&alphaf,&betaf);   
    }
    fclose(fp);
    *alpha = (SSAL_real_t)alphaf;
    *beta = (SSAL_real_t)betaf;
    return 0;
}

void GetCmdLineArgs(int argc,char **argv)
{
    int i;
    for (i=1;i<argc;i++)
    {
        if (!strcmp("-e",argv[i])) /*target error*/
        {
            ml_p.target_RMSE = (SSAL_real_t) atof(argv[++i])/sqrt(4);           
        }
        else if (!strcmp("-f",argv[i]))
        {
            model_filename = argv[++i];   
            data_filename = argv[++i];   
        }
        else if (!strcmp("-M",argv[i]))
        {
            ml_p.K = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-eps0",argv[i]))
        {
            ml_p.eps0 = (SSAL_real_t)atof(argv[++i]);
        }
        else if (!strcmp("-div",argv[i]))
        {
           D = (int)atoi(argv[++i]);
        }
        else if (!strcmp("-param",argv[i]))
        {
            ml_p.alpha = (SSAL_real_t) atof(argv[++i]);
            ml_p.beta = (SSAL_real_t) atof(argv[++i]);
            ml_p.L = ceil((log(1.0/ml_p.target_RMSE)/log(ml_p.K))/ml_p.alpha);
        }
        else if (!strcmp("-paramf",argv[i]))
        {
            char *parafile;
            unsigned int d;
            parafile = argv[++i];
            d = (unsigned int)atoi(argv[++i]);
            getCovergenceParams(parafile,d,&(ml_p.alpha), &(ml_p.beta));
            ml_p.L = ceil((log(1.0/ml_p.target_RMSE)/log(ml_p.K))/ml_p.alpha);
        }
        else if (!strcmp("-ml",argv[i]))
        {
            ml = 1;
        }
        else if (!strcmp("-o",argv[i]))
        {
            outputfilename = argv[++i];
            summaryfilename = argv[++i];
        }
        else if (!strcmp("-h",argv[i]))
        {
            fprintf(stderr,"%s: -e [target error] -f [modelfile datafile] -M [scale factor] -eps0 [base epsilon] -div [num indicator functions per axis] -param [alpha beta] -o [output summary] [-ml] -s [smoothfactor scalefactor] \n",argv[0]);
        }
        else if (!strcmp("-s",argv[i]))
        {
            smoothfactor = (SSAL_real_t)atof(argv[++i]);
            scalefactor = (SSAL_real_t)atof(argv[++i]);
        }
        else if (!strcmp("-L",argv[i]))
        {
            ml_p.L = atoi(argv[++i]);
        }
        else if (!strcmp("-a",argv[i]))
        {
            /*use theoretical bounds*/
            ml_p.presample = 0;
        }
    }
}

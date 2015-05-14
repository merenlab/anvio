#ifndef NMGS_H
#define NMGS_H

typedef struct s_Params
{
    /*seed*/
    unsigned long int lSeed;
    /*min change VBL*/
    double dEpsilon;
    /*maximum no. iterations*/
    int nMaxIter;
    /*initial cluster size*/
    int nKStart;
} t_Params;


typedef struct s_Data
{
    int nN;

    int nD;

    double **aadX;
} t_Data;

typedef struct s_VBParams
{
    /*scale for mean prior*/
    double dBeta0;

    /*Wishart degrees of freedom*/
    double dNu0;

    /*Inverse! of the Wishart scale precision-matrix*/
    gsl_matrix *ptInvW0;

    /*Log Wishart normalisation*/
    double dLogWishartB;

} t_VBParams;


typedef struct s_Cluster
{
    /*output file for convergence if not null*/
    char *szCOutFile;
    /*parameters for variational Bayes*/
    t_VBParams *ptVBParams;
    /*start seed*/
    unsigned long lSeed;
    /* maximum no. iterations*/
    int nMaxIter;
    /* min. change in VBL bound*/
    double dEpsilon;
    /*thread index*/
    int nThread;
    /*pointer to data*/
    t_Data *ptData;
    /*number of data points*/
    int nN;
    /*size no. of clusters allocated*/
    int nKSize;
    /*number of clusters*/
    int nK;
    /*number of dimensions*/
    int nD;
    /*variational lower bound*/
    double dVBL;
    /*Means*/
    double **aadMu;
    /*Scaled weight Bishop 10.60*/
    double *adBeta;
    /*Scaled means Bishop 10.61*/
    double **aadM;
    /*sample covariance matrix for each cluster storing this helps with lower bound calcn*/
    gsl_matrix **aptCovar;
    /*Inverse regularised variances Bishop 10.62*/
    gsl_matrix **aptSigma;
    /*Bishop 10.63*/
    double *adNu;
    /*Responsibilities*/
    double **aadZ;
    /*log-Matrix determinants*/
    double *adLDet;
    /*mixture weights*/
    double *adPi;
    /*assigned cluster for each data point*/
    int *anMaxZ;
    /*frequencies for each cluster*/
    int *anW;
} t_Cluster;


#define DELIM ",\n"
#define MAX_LINE_LENGTH   1048576
#define MAX_FILE_NAME_LENGTH 1024
#define MAX_WORD_LENGTH   128
#define DEF_FILE_STUB     "debug_out"

#define TRUE  1
#define FALSE 0

#define NOT_SET -1

/*Default parameters*/
#define DEF_BETA0        1.0e-3

#define MIN_Z            1.0e-6
#define MIN_PI           0.1 /*Unormalised*/
#define MIN_COVAR        0.001

#define N_RTHREADS       2
#define R_PRIME          1009


#define DEF_EPSILON      1.0e-6
#define DEF_MAX_ITER     500
#define DEF_SEED         1l

/*user defines*/
int driver(double *adX, int nN, int nD, int *anAssign, int nKStart, unsigned long lSeed, int nMaxIter, double dEpsilon, int debug);

void generateInputData(double *adX, int nN, int nD, t_Data *ptData);

void destroyData(t_Data *ptData);

void calcSampleVar(t_Data *ptData,double *adVar, double *adMu);

void setVBParams(t_VBParams *ptVBParams, t_Data *ptData);

void* fitEM(void *pvCluster);

void* runRThreads(void *pvpDCluster);

void allocateCluster(t_Cluster *ptCluster, int nN, int nK, int nD, t_Data *ptData, long lSeed, int nMaxIter, double dEpsilon, char *szCOutFile);

void destroyCluster(t_Cluster* ptCluster);

void compressCluster(t_Cluster *ptCluster);

double decomposeMatrix(gsl_matrix *ptSigmaMatrix, int nD);

void performMStep(t_Cluster *ptCluster, t_Data *ptData);

void initKMeans(gsl_rng *ptGSLRNG, t_Cluster *ptCluster, t_Data *ptData);

double calcVBL(t_Cluster* ptCluster);

void gmmTrainVB(t_Cluster *ptCluster, t_Data *ptData);

double dLogWishartB(gsl_matrix *ptInvW, int nD, double dNu, int bInv);

void updateMeans(t_Cluster *ptCluster, t_Data *ptData);

double dWishartExpectLogDet(gsl_matrix *ptW, double dNu, int nD);

void calcZ(t_Cluster* ptCluster, t_Data *ptData);

void calcCovarMatrices(t_Cluster *ptCluster, t_Data *ptData);

double calcDist(double* adX, double *adMu, int nD);

#endif

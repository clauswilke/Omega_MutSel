/* SJS 9/9/14 - 9/10/14.
Hyphy inference for an "experimental" dataset. Name of file indicates the mutation scheme.
Perform 7 inferences for each of the following parameterizations: Fequal, F61_true, F61_data, F3x4_true, F61_data, Fnuc_true, Fnuc_data. The _data refers to empirical frequencies, whereas _true refers to frequencies in absence of selection. 
Also note that Fnuc is not so much a frequency parameterization, but actually a "new" model.
*/


global w;
global k;
global t;

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
#include "GY94.mdl"; // Basic GY94 rate matrix
#include "fnuc.mdl"; // Custom Fnuc matrices for this run

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


/* Set up frequencies. MANY OF THESE WERE HARD-CODED IN WHEN THIS FILE WAS CREATED!!!*/

Fequal = {{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623}};

F61_true = INSERT_F61_TRUE

F61_data = INSERT_F61_DATA

F3x4_true = INSERT_F3X4_TRUE

F3x4_data = INSERT_F3X4_DATA


/* Optimize likelihoods for each frequency specification */

////////////// FEQUAL FREQUENCIES //////////////
Model MyModel = (GY94, Fequal, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn1 = (filt_data, Tree01);
Optimize (paramValues, LikFn1);
fprintf ("fequal_hyout.txt", LikFn1);


////////////// F61_TRUE FREQUENCIES //////////////
Model MyModel = (GY94, F61_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn2 = (filt_data, Tree01);
Optimize (paramValues, LikFn2);
fprintf ("f61_true_hyout.txt", LikFn2);



////////////// F61_DATA FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, F61_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn3 = (filt_data, Tree01);
Optimize (paramValues, LikFn3);
fprintf ("f61_data_hyout.txt", LikFn3);



////////////// F3x4_TRUE FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, F3x4_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn4 = (filt_data, Tree01);
Optimize (paramValues, LikFn4);
fprintf ("f3x4_true_hyout.txt", LikFn4);


////////////// F3x4_DATA FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, F3x4_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn5 = (filt_data, Tree01);
Optimize (paramValues, LikFn5);
fprintf ("f3x4_data_hyout.txt", LikFn5);



Fones =  {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
////////////// Fnuc_TRUE MODEL //////////////
global w;
global k;
global t;
Model MyModel = (GY94_Fnuc_true, Fones, 0); // Using 0 as last argument means that the matrix will *not* be multipled by frequencies, but just in case it is, we provide Fones (all entries are 1, so multiplication is basically..not)
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn6 = (filt_data, Tree01);
Optimize (paramValues, LikFn6);
fprintf ("fnuc_true_hyout.txt", LikFn6);


////////////// Fnuc_DATA MODEL //////////////
global w;
global k;
global t;
Model MyModel = (GY94_Fnuc_data, Fones, 0); // Using 0 as last argument means that the matrix will *not* be multipled by frequencies, but just in case it is, we provide Fones (all entries are 1, so multiplication is basically..not)
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn7 = (filt_data, Tree01);
Optimize (paramValues, LikFn7);
fprintf ("fnuc_data_hyout.txt", LikFn7);

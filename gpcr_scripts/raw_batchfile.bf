/* SJS. 
Hyphy inference for GPCR data sets
*/

global w; global k; global t;

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;
#include "GY94.mdl"; // Basic GY94 rate matrix
#include "fnuc.mdl"; // Custom Fnuc matrices for this run

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


/* Set up frequencies. Note that these were all hard-coded in when the file was created via the script Omega_MutSel/np_scripts/prefs_to_freqs.py */

F61 = INSERT_F61

F1x4 = INSERT_F1X4

F3x4 = INSERT_F3X4

CF3x4 = INSERT_CF3X4



/* Optimize likelihoods for each frequency specification */
////////////// F61 FREQUENCIES //////////////
Model MyModel = (GY94, F61, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn1 = (filt_data, Tree01);
Optimize (paramValues, LikFn1);
fprintf ("f61_hyout.txt", LikFn1);


////////////// F1x4 FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F1x4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn3 = (filt_data, Tree01);
Optimize (paramValues, LikFn3);
fprintf ("f1x4_hyout.txt", LikFn3);


////////////// F3x4_TRUE FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F3x4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn5 = (filt_data, Tree01);
Optimize (paramValues, LikFn5);
fprintf ("f3x4_hyout.txt", LikFn5);


////////////// CF3x4 FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, CF3x4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn8 = (filt_data, Tree01);
Optimize (paramValues, LikFn8);
fprintf ("cf3x4_hyout.txt", LikFn8);


////////////// Fnuc MODEL //////////////
Fones =  {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
global w; global k; global t;
Model MyModel = (Fnuc, Fones, 0); // Using 0 as last argument means that the matrix will *not* be multipled by frequencies, but just in case it is, we provide Fones (all entries are 1, so multiplication is basically..not)
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn9 = (filt_data, Tree01);
Optimize (paramValues, LikFn9);
fprintf ("fnuc_hyout.txt", LikFn9);

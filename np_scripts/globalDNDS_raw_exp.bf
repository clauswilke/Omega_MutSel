/* SJS. 
Hyphy inference for an "experimental" dataset. Name of file indicates the mutation scheme.
Perform 10 total inferences, one for each of the following parameterizations: F61_true, F61_data, F1x4_true, F1x4_data, F3x4_true, F3x4_data, CF3x4_true, CF3x4_data, Fnuc_true, Fnuc_data. The _data refers to empirical frequencies, whereas _true refers to frequencies in absence of selection. 
Also note that Fnuc is not so much a frequency parameterization, but actually a "new"(ish? it's actually what should have been the original) model.
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

F61_true = INSERT_F61_TRUE

F61_data = INSERT_F61_DATA

F1x4_true = INSERT_F1X4_TRUE

F1x4_data = INSERT_F1X4_DATA

F3x4_true = INSERT_F3X4_TRUE

F3x4_data = INSERT_F3X4_DATA

CF3x4_true = INSERT_CF3X4_TRUE

CF3x4_data = INSERT_CF3X4_DATA

/* Optimize likelihoods for each frequency specification */


////////////// F61_TRUE FREQUENCIES //////////////
Model MyModel = (GY94, F61_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn1 = (filt_data, Tree01);
Optimize (paramValues, LikFn1);
fprintf ("f61_true_hyout.txt", LikFn1);



////////////// F61_DATA FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F61_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn2 = (filt_data, Tree01);
Optimize (paramValues, LikFn2);
fprintf ("f61_data_hyout.txt", LikFn2);


////////////// F1x4_TRUE FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F1x4_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn3 = (filt_data, Tree01);
Optimize (paramValues, LikFn3);
fprintf ("f1x4_true_hyout.txt", LikFn3);


////////////// F1x4_DATA FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F1x4_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn4 = (filt_data, Tree01);
Optimize (paramValues, LikFn4);
fprintf ("f1x4_data_hyout.txt", LikFn4);


////////////// F3x4_TRUE FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F3x4_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn5 = (filt_data, Tree01);
Optimize (paramValues, LikFn5);
fprintf ("f3x4_true_hyout.txt", LikFn5);


////////////// F3x4_DATA FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, F3x4_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn6 = (filt_data, Tree01);
Optimize (paramValues, LikFn6);
fprintf ("f3x4_data_hyout.txt", LikFn6);

////////////// CF3x4_TRUE FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, CF3x4_true, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn7 = (filt_data, Tree01);
Optimize (paramValues, LikFn7);
fprintf ("cf3x4_true_hyout.txt", LikFn7);


////////////// CF3x4_DATA FREQUENCIES //////////////
global w; global k; global t;
Model MyModel = (GY94, CF3x4_data, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn8 = (filt_data, Tree01);
Optimize (paramValues, LikFn8);
fprintf ("cf3x4_data_hyout.txt", LikFn8);


Fones =  {{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1},{1}};
////////////// Fnuc_TRUE MODEL //////////////
global w; global k; global t;
Model MyModel = (Fnuc_true, Fones, 0); // Using 0 as last argument means that the matrix will *not* be multipled by frequencies, but just in case it is, we provide Fones (all entries are 1, so multiplication is basically..not)
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn9 = (filt_data, Tree01);
Optimize (paramValues, LikFn9);
fprintf ("fnuc_true_hyout.txt", LikFn9);


////////////// Fnuc_DATA MODEL //////////////
global w; global k; global t;
Model MyModel = (Fnuc_data, Fones, 0); // Using 0 as last argument means that the matrix will *not* be multipled by frequencies, but just in case it is, we provide Fones (all entries are 1, so multiplication is basically..not)
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn10 = (filt_data, Tree01);
Optimize (paramValues, LikFn10);
fprintf ("fnuc_data_hyout.txt", LikFn10);


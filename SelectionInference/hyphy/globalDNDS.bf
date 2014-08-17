/* Hyphy inference across 4 sets of equilibrium freqs: equal, F61, F3x4, CF3x4. */

global w;
global k;
global t;

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;

/* Include relevant functions */
#include "matrices.mdl"; //rate matrices
#include "GY94_Header.ibf";

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");



/* Set up frequencies. 1/61, F61, F3x4, CF3x4 */

equal_CodonFreqs = {{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623},{0.016393442623}};

F61_CodonFreqs = DATAFREQS

HarvestFrequencies(posFreqs,filt_data,3,1,1);
F3x4_CodonFreqs = BuildCodonFrequencies(posFreqs);

posFreqs_cf3x4 = CF3x4(posFreqs, "TAA,TAG,TGA");
CF3x4_CodonFreqs = BuildCodonFrequencies(posFreqs_cf3x4);


/* Optimize likelihoods for each frequency specification */

////////////// EQUAL FREQUENCIES //////////////
Model MyModel = (GY94, equal_CodonFreqs, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn1 = (filt_data, Tree01);
Optimize (paramValues, LikFn1);
fprintf ("equal_hyout.txt", LikFn1);


////////////// DATA FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, F61_CodonFreqs, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn2 = (filt_data, Tree01);
Optimize (paramValues, LikFn2);
fprintf ("f61_hyout.txt", LikFn2);

////////////// F3x4 FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, F3x4_CodonFreqs, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn3 = (filt_data, Tree01);
Optimize (paramValues, LikFn3);
fprintf ("f3x4_hyout.txt", LikFn3);

////////////// CF3x4 FREQUENCIES //////////////
global w;
global k;
global t;
Model MyModel = (GY94, CF3x4_CodonFreqs, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn4 = (filt_data, Tree01);
Optimize (paramValues, LikFn4);
fprintf ("cf3x4_hyout.txt", LikFn4);


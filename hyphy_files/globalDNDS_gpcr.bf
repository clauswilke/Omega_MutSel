/* SJS. 
Hyphy inference for GPCR data sets
*/

global w; global k;

LIKELIHOOD_FUNCTION_OUTPUT = 1;
RANDOM_STARTING_PERTURBATIONS = 1;
OPTIMIZATION_PRECSION = 0.00000001;
#include "CF3x4.bf"; // to make the CF3x4 frequencies
#include "GY94.mdl"; // Basic GY94 rate matrix
#include "fnuc.mdl"; // Custom Fnuc matrices for this run

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");

HarvestFrequencies(freqVector,filt_data,1,1,1);
HarvestFrequencies(freqVector_pos,filt_data,3,1,1);
HarvestFrequencies(freqVector_61,filt_data,3,3,1);
F61 = Transpose(freqVector_61[Transpose(_Genetic_Code["_MATRIX_ELEMENT_VALUE_!=10"])]); /* cough, cough, cough....HACK! */ // <- sergei, not me! :)
F1x4 = BuildCodonFrequencies(freqVector);
F3x4 = BuildCodonFrequencies(freqVector_pos);
CF3x4_ = BuildCodonFrequencies(CF3x4(freqVector_pos, "TAA,TAG,TGA"));


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
global w; global k;
Model MyModel = (GY94, F1x4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn3 = (filt_data, Tree01);
Optimize (paramValues, LikFn3);
fprintf ("f1x4_hyout.txt", LikFn3);


////////////// F3x4_TRUE FREQUENCIES //////////////
global w; global k;
Model MyModel = (GY94, F3x4, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn5 = (filt_data, Tree01);
Optimize (paramValues, LikFn5);
fprintf ("f3x4_hyout.txt", LikFn5);


////////////// CF3x4 FREQUENCIES //////////////
global w; global k;
Model MyModel = (GY94, CF3x4_, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn8 = (filt_data, Tree01);
Optimize (paramValues, LikFn8);
fprintf ("cf3x4_hyout.txt", LikFn8);


////////////// Fnuc_pos MODEL //////////////
global w; global k;
Model MyModel = (Fnuc_pos, F3x4, 0);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn9 = (filt_data, Tree01);
Optimize (paramValues, LikFn9);
fprintf ("fnuc_pos_hyout.txt", LikFn9);

////////////// Fnuc_glob MODEL //////////////
global w; global k;
Model MyModel = (Fnuc_glob, F1x4, 0);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn9 = (filt_data, Tree01);
Optimize (paramValues, LikFn9);
fprintf ("fnuc_glob_hyout.txt", LikFn9);
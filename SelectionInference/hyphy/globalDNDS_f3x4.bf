global a;
global b;
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

HarvestFrequencies(posFreqs,filt_data,3,1,1);
codonFreq = BuildCodonFrequencies(posFreqs);


Model MyModel = (GY94, codonFreq, 1);
UseModel (USE_NO_MODEL);
UseModel(MyModel);
Tree    Tree01 = DATAFILE_TREE;
LikelihoodFunction  LikFn = (filt_data, Tree01);
Optimize (paramValues, LikFn);
fprintf (stdout, LikFn);

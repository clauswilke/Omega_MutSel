#include "../SelectionInference/hyphy/GY94_Header.ibf";


OPTIMIZATION_PRECISION = 0.00001;

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");

HarvestFrequencies(posFreqs,filt_data,3,1,1);

f3x4 = BuildCodonFrequencies(posFreqs);
fprintf("f3x4.txt", f3x4);

posFreqs_cf3x4 = CF3x4(posFreqs, "TAA,TAG,TGA"); // converts f3x4 positional nucleotide frequencies to cf3x4 positional nucleotide frequencies
cf3x4 = BuildCodonFrequencies(posFreqs_cf3x4);
fprintf("cf3x4.txt", cf3x4);
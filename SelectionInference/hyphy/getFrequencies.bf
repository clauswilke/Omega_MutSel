/* Just grab F3x4, CF3x4 frequencies */


/* Include relevant functions */
#include "GY94_Header.ibf";

/* Read in the data */
DataSet	raw_data = ReadDataFile("temp.fasta");

/* Filter the data to find and remove any stop codons*/
DataSetFilter   filt_data = CreateFilter(raw_data,3,"", "","TAA,TAG,TGA");


HarvestFrequencies(posFreqs,filt_data,3,1,1);
F3x4_CodonFreqs = BuildCodonFrequencies(posFreqs);
fprintf("f3x4.txt", F3x4_CodonFreqs);

posFreqs_cf3x4 = CF3x4(posFreqs, "TAA,TAG,TGA");
CF3x4_CodonFreqs = BuildCodonFrequencies(posFreqs_cf3x4);
fprintf("cf3x4.txt", CF3x4_CodonFreqs);

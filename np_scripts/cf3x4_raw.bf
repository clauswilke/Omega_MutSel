#include "/usr/local/lib/hyphy/TemplateBatchFiles/TemplateModels/CF3x4.bf"; 
VERBOSITY_LEVEL=0;
OPTIMIZATION_PRECISION = 0.000000001;
pos_freqs = INSERT_POS_FREQS
fprintf (stdout, CF3x4(pos_freqs, "TAA,TAG,TGA"));

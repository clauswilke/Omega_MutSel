Omega_MutSel
============

Repository for "The relationship between dN/dS and scaled selection coefficients", Stephanie J. Spielman and Claus O. Wilke.
All code written by SJS (contact at stephanie.spielman@gmail.com).

Contents: 

*datasets/  Contains tab-delimited summary files for simulated datasets.  All simulated alignments available upon request from stephanie.spielman@gmail.com.
 *nobias.txt -> simulations with symmetric mutation rates and in which synonymous codons all have same fitness (Figure 1A, Figure 2B)
 *bias.txt   -> simulations with symmetric mutation rates and in synonymous codons have different fitness (Figure 1B, Figure 2B)
 *conv.txt   -> simulations to demonstrate convergence of omega to dN/dS (Figure 2C)
 *np.txt     -> simulations which use experimental NP amino acid fitness data (Bloom 2014) in combination with NP mutation rates (Bloom 2014) (Figure 3, Tables 1,S1,S2, and Figure S1)
 *yeast.txt  -> simulations which use experimental NP amino acid fitness data (Bloom 2014) in combination with yeast mutation rates (Zhu 2014) (Figure 3, Tables 1,S1,S2, and Figure S1)
 *polio.txt  -> simulations which use experimental NP amino acid fitness data (Bloom 2014) in combination with polio virus mutation rates (Acevedo 2014) (Figure 3, Tables 1,S1,S2, and Figure S1)

*scripts/ Contains scripts used in analysis. [NOTE: all simulated alignments were created using a custom sequence simulation library, available from https://github.com/sjspielman/pyvolve. See within for details.]

 *experimental_data/
  *nucleoprotein_amino_preferences.txt -> This file corresponds exactly to Table S1 of Bloom 2014. Gives amino acid preference/fitness data for each of the 498 positions in NP. Each row is a position, and each column is the amino acid preference (alphabetical)
  *np_codon_eqfreqs.txt                -> Contains codon equilibrium frequencies computed from NP preference data and NP mutation rates. Each row is a position, and values are alphabetical (first column is AAA, second column is AAC, etc.). Generated by scripts/np_scripts/prefs_to_freqs.py .
  *yeast_codon_eqfreqs.txt             -> Contains codon equilibrium frequencies computed from yeast preference data and yeast mutation rates. Each row is a position, and values are alphabetical (first column is AAA, second column is AAC, etc.). Generated by scripts/yeast_scripts/prefs_to_freqs.py .
  *polio_codon_eqfreqs.txt             -> Contains codon equilibrium frequencies computed from polio preference data and polio mutation rates. Each row is a position, and values are alphabetical (first column is AAA, second column is AAC, etc.). Generated by scripts/polio_scripts/prefs_to_freqs.py .

 *np_scripts/
  *prefs_to_freqs.py      -> Compute equilibrium codon frequencies for a variety of frequency parameterizations from experimental NP amino acid fitness data in combination with either NP, yeast, or polio mutation rates. See script for full description.
  *globalDNDS_raw_exp.bf  -> Template batchfile used by prefs_to_freqs.py to create hyphy_files/globalDNDS_{np/yeast/polio}.bf

 *simulation_scripts/   Scripts in this directory were created to run specifically on The University of Texas at Austin's Center for Computational Biology and Bioinformatics cluster, "Phylocluster". All files w/ extension ".qsub" are job submission scripts corresponding to a particular python script, such that <xyz>.qsub goes with run_<xyz>.py.
  *
  *sim_nyp.py            -> simulate alignments which use NP amino acid fitness data and either NP, yeast, or polio mutation rates
  *run_nyp.py             -> infer dN/dS, omega for NP, yeast, or polio datasets
  *run_siminf.py          -> simulate alignments, infer dN/dS and omega for the "codon bias" and "no codon bias" sets
  *run_convergence.py     -> simulate alignmets, infer dN/dS and omega to demonstrate omega convergence with data sets of increasing size
  *functions_simandinf.py -> contains functions used by scripts in this directory.

 *hyphy_files/        Contains files used in HYPHY inference.
  *globalDNDS_fequal.bf -> hyphy batchfile to infer omega according to GY94 M0 model with Fequal (1/61 for all codons) frequency parameterization. Used to determined omega for nobias.txt, bias.txt, conv.txt .
  *globalDNDS_np.bf     -> hyphy batchfile to infer omega for simulations with experimental NP amino acid fitness data and NP mutation rates, according to a variety of frequency parameterizations. Created by scripts/np_scripts/prefs_to_freqs.py .
  *globalDNDS_yeast.bf  -> hyphy batchfile to infer omega for simulations with experimental NP amino acid fitness data and yeast mutation rates, according to a variety of frequency parameterizations. Created by scripts/np_scripts/prefs_to_freqs.py .
  *globalDNDS_polio.bf  -> hyphy batchfile to infer omega for simulations with experimental NP amino acid fitness data and polio mutation rates, according to a variety of frequency parameterizations. Created by scripts/np_scripts/prefs_to_freqs.py .
  *CF3x4.bf             -> hyphy batchfile used in conjunction with globalDNDS_{np,polio,yeast}.bf to compute CF3x4 equilibrium codon frequencies.
  *GY94.mdl             -> contains GY94 rate matrix
  *MG_np.mdl            -> contains MG1 and MG3 matrices for NP mutation rates 
  *MG_yeast.mdl         -> contains MG1 and MG3 matrices for yeast mutation rates 
  *MG_polio.mdl         -> contains MG1 and MG3 matrices for polio mutation rates 












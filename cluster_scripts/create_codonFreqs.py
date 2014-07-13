

import sys
sys.path.append('/home/sjs3495/MutSel/Simulator/src')
from functions_simandinf import *


for i in range(1,51):
    beta = rn.uniform(0.5,3.5)
    f, num_pref_aa, gc_content = setFreqs('codonFreqs/codonFreq'+str(i)+'.txt', beta, 0.0, 1.0) # last 2 args are gc min, gc max




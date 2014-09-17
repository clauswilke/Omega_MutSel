_Genetic_Code = {
                     {14,13,14,13,7,7,7,7,19, 5,19, 5,2,2,3,2,
                      12,11,12,11,6,6,6,6,19,19,19,19,1,1,1,1,
                      16,15,16,15,8,8,8,8,20,20,20,20,4,4,4,4,
                      10,9, 10,9, 5,5,5,5,10,17,18,17,1,0,1,0}
                };
Codon_Index = {
				{14,13,14,13, 7, 7, 7, 7,19, 5,19, 5, 2, 2, 3, 2,
				 12,11,12,11, 6, 6, 6, 6,19,19,19,19, 1, 1, 1, 1,
				 16,15,16,15, 8, 8, 8, 8,20,20,20,20, 4, 4, 4, 4,
				     9,    9, 5, 5, 5, 5,   17,18,17, 1, 0, 1, 0}
			}; 

Codon_Index_ENG = {

				{"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
				  "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
				 "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
				 "TAC",      "TAT","TCA","TCC","TCG","TCT",      "TGC","TGG","TGT","TTA","TTC","TTG","TTT"}
           	   };

function BuildCodonFrequencies (obsF)
{
        PIStop = 1.0;
        result = {61,1};
        hshift = 0;

        for (h=0; h<64; h=h+1)
        {
                first = h$16;
                second = h%16$4;
                third = h%4;
                if (_Genetic_Code[h]==10) 
                {
                        hshift = hshift+1;
                        PIStop = PIStop-obsF[first][0]*obsF[second][1]*obsF[third][2];
                        continue; 
                }
                result[h-hshift][0]=obsF[first][0]*obsF[second][1]*obsF[third][2];
                
        }

        return result*(1.0/(PIStop))
        
}

function CF3x4 (observed3x4,stopCodons)
{
	global p11 :< 1; p11:>0;
	global p12 :< 1; p12:>0;
	global p13 :< 1; p13:>0;
	global p21 :< 1; p21:>0;
	global p22 :< 1; p22:>0;
	global p23 :< 1; p23:>0;
	global p31 :< 1; p31:>0;
	global p32 :< 1; p32:>0;
	global p33 :< 1; p33:>0;
	
	N = {4,3};
	
	n1A : = p11;
	p11 = observed3x4[0][0];
	n1C : = (1-p11)*p12;
	p12 = observed3x4[1][0]/(1-p11);
	n1G : = (1-p11)*(1-p12)*p13;
	p13 = observed3x4[2][0]/(1-p11)/(1-p21);
	n1T : = (1-p11)*(1-p12)*(1-p13);

	n2A : = p21;
	p21 = observed3x4[0][1];
	n2C : = (1-p21)*p22;
	p22 = observed3x4[1][1]/(1-p21);
	n2G : = (1-p21)*(1-p22)*p23;
	p23 = observed3x4[2][1]/(1-p21)/(1-p22);
	n2T : = (1-p21)*(1-p22)*(1-p23);

	n3A : = p31;
	p31 = observed3x4[0][2];
	n3C : = (1-p31)*p32;
	p32 = observed3x4[1][2]/(1-p31);
	n3G : = (1-p31)*(1-p32)*p33;
	p33 = observed3x4[2][2]/(1-p31)/(1-p32);
	n3T : = (1-p31)*(1-p32)*(1-p33);
	
	stopCount 		= (Abs(stopCodons)+1)$4;
	charMap   		= {"A":0,"C":1,"G":2, "T": 3};
	revMap			= {{"A","C","G","T"}};
	stopComposition = {stopCount,3};
	
	SDef =""; SDef * 128; SDef * "1";
	
	for (i = 0; i < stopCount; i = i+1)
	{
		SDef * "-";
		for (j = 0; j < 3; j = j+1)
		{
			stopComposition[i][j] = charMap[stopCodons[4*i+j]];
			if (j)
			{
				SDef * "*";
			}
			SDef * ("n"+(j+1)+stopCodons[4*i+j]);
		}
	}
	
	SDef * 0;
	
	ExecuteCommands ("global S:=`SDef`");
	
	SDef = {4,3};
	for (i = 0; i<4; i=i+1)
	{
		for (j = 0; j<3; j=j+1)
		{
			SDef[i][j] = "";
		}
	}
	
	for (k = 0; k < stopCount; k = k+1)
	{
		SDef[stopComposition[k][0]][0] = SDef[stopComposition[k][0]][0] + 
										 "-" +
										 "n2" + revMap[stopComposition[k][1]] +
										 "*n3" + revMap[stopComposition[k][2]];
		SDef[stopComposition[k][1]][1] = SDef[stopComposition[k][1]][1] + 
										 "-" +
										 "n1" + revMap[stopComposition[k][0]] +
										 "*n3" + revMap[stopComposition[k][2]];
		SDef[stopComposition[k][2]][2] = SDef[stopComposition[k][2]][2] + 
										 "-" +
										 "n1" + revMap[stopComposition[k][0]] +
										 "*n2" + revMap[stopComposition[k][1]];
	}

	for (i = 0; i<4; i=i+1)
	{
		for (j = 0; j<3; j=j+1)
		{
			if (Abs (SDef[i][j]))
			{
				ExecuteCommands ("N[i][j] := n"+(j+1)+revMap[i] + "*(1" + SDef[i][j] + ")/S;");
			}
			else
			{
				ExecuteCommands ("N[i][j] := n"+(j+1)+revMap[i] + "/S;");			
			}
		}
	}

	calc = 0;
	
	norm1 = {3,1}["1"];
	norm2 = {1,4}["1"];
	
	
	Optimize (res, _CF3x4_minimizer(p11,p12,p13,p21,p22,p23,p31,p32,p33));
	
	return {{n1A__,n2A__,n3A__}{n1C__,n2C__,n3C__}{n1G__,n2G__,n3G__}{n1T__,n2T__,n3T__}};
}

function _CF3x4_minimizer (p11,p12,p13,p21,p22,p23,p31,p32,p33)
{
	calc = calc + 1;
	error = N-observed3x4;
	error = error$error;
	return  -(norm2*error*norm1)[0];
}

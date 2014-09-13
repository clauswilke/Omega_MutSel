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

F61_true = {{0.0415420976134},{0.019859659982},{0.0199226668741},{0.0418410026486},{0.019859659982},{0.00952086207663},{0.00952773729476},{0.0199841786885},{0.0199226668741},{0.00952773729476},{0.00960027733595},{0.0200172808268},{0.0418410026486},{0.0199841786885},{0.0200172808268},{0.0418563181783},{0.0212103555792},{0.00960702697888},{0.0100649760394},{0.0202818863501},{0.00960702697888},{0.00456888973382},{0.00459672580976},{0.00960216068013},{0.0100649760394},{0.00459672580976},{0.00472709614948},{0.00965950242228},{0.0202818863501},{0.00960216068013},{0.00965950242228},{0.0201068649502},{0.0189747676791},{0.00936938936388},{0.00917571428658},{0.0199620548488},{0.00936938936388},{0.00452903599734},{0.00451834486019},{0.0095442790625},{0.00917571428658},{0.00451834486019},{0.0045759964321},{0.00956040435674},{0.0199620548488},{0.0095442790625},{0.00956040435674},{0.0200222354502},{0.0189970705278},{0.0428338830373},{0.0189970705278},{0.00937204102203},{0.00925761483226},{0.0200028931788},{0.00925761483226},{0.0105071511614},{0.0202173800273},{0.0428338830373},{0.0200028931788},{0.0202173800273},{0.0420783446859}};

F61_data = {{0.0425878663474},{0.0199078643132},{0.0206080077075},{0.0410957486733},{0.0211288327264},{0.0101737510099},{0.0102015024614},{0.0211244049354},{0.0281069065528},{0.0113126700256},{0.0135933116374},{0.0233415011914},{0.038607136559},{0.0183427368963},{0.0208081941436},{0.0378467556072},{0.0212454327197},{0.00799840204442},{0.0102324620869},{0.0165134584925},{0.00727186274848},{0.00347178512148},{0.00349349462},{0.00721584392773},{0.0144718989509},{0.00684780307396},{0.00682240822034},{0.0139887873167},{0.0205386820914},{0.00972918740468},{0.0097893064169},{0.0198883449017},{0.0190032600421},{0.00794583686669},{0.00925901840985},{0.0165008466301},{0.0130669480887},{0.00632403196468},{0.00631782177592},{0.0132104311994},{0.012572795276},{0.00617384774217},{0.00622587915583},{0.0128712793025},{0.0193448210618},{0.00925162260606},{0.00926846715473},{0.0191727272757},{0.0136570776602},{0.03037549023},{0.0221599938113},{0.010960463496},{0.0108055448481},{0.0230642639634},{0.00804200517678},{0.00718219624967},{0.017182851242},{0.0429633607131},{0.0156993564494},{0.0205206463864},{0.0325707642958}};

F1x4_true = {{0.0366337204149},{0.0197450646254},{0.0182477196333},{0.0384277729969},{0.0197450646254},{0.0106423145846},{0.00983526650707},{0.0207120339578},{0.0182477196333},{0.00983526650707},{0.00908942002191},{0.019141359923},{0.0384277729969},{0.0207120339578},{0.019141359923},{0.0403096851965},{0.0197450646254},{0.0106423145846},{0.00983526650707},{0.0207120339578},{0.0106423145846},{0.00573605920604},{0.00530107154257},{0.0111634975751},{0.00983526650707},{0.00530107154257},{0.00489907068426},{0.0103169261658},{0.0207120339578},{0.0111634975751},{0.0103169261658},{0.0217263583993},{0.0182477196333},{0.00983526650707},{0.00908942002191},{0.019141359923},{0.00983526650707},{0.00530107154257},{0.00489907068426},{0.0103169261658},{0.00908942002191},{0.00489907068426},{0.00452755511345},{0.00953455355665},{0.019141359923},{0.0103169261658},{0.00953455355665},{0.0200787641998},{0.0207120339578},{0.0403096851965},{0.0207120339578},{0.0111634975751},{0.0103169261658},{0.0217263583993},{0.0103169261658},{0.00953455355665},{0.0200787641998},{0.0403096851965},{0.0217263583993},{0.0200787641998},{0.0422837597373}};

F1x4_data = {{0.037367588363},{0.0204382258569},{0.0206300848715},{0.0360963515048},{0.0204382258569},{0.011178700432},{0.0112836378402},{0.0197429220611},{0.0206300848715},{0.0112836378402},{0.0113895603235},{0.0199282540757},{0.0360963515048},{0.0197429220611},{0.0199282540757},{0.0348683618353},{0.0204382258569},{0.011178700432},{0.0112836378402},{0.0197429220611},{0.011178700432},{0.006114197202},{0.00617159278315},{0.0107984035855},{0.0112836378402},{0.00617159278315},{0.00622952715175},{0.0108997710469},{0.0197429220611},{0.0107984035855},{0.0108997710469},{0.0190712723424},{0.0206300848715},{0.0112836378402},{0.0113895603235},{0.0199282540757},{0.0112836378402},{0.00617159278315},{0.00622952715175},{0.0108997710469},{0.0113895603235},{0.00622952715175},{0.00628800536554},{0.0110020900715},{0.0199282540757},{0.0108997710469},{0.0110020900715},{0.019250299404},{0.0197429220611},{0.0348683618353},{0.0197429220611},{0.0107984035855},{0.0108997710469},{0.0190712723424},{0.0108997710469},{0.0110020900715},{0.019250299404},{0.0348683618353},{0.0190712723424},{0.019250299404},{0.0336821480952}};

F3x4_true = {{0.0355258267054},{0.0202241751754},{0.018243496951},{0.0430053423838},{0.0202241751754},{0.0115132369731},{0.0103856746589},{0.0244821207135},{0.018243496951},{0.0103856746589},{0.00936854147721},{0.0220844356181},{0.0430053423838},{0.0244821207135},{0.0220844356181},{0.0520595759496},{0.0173564056366},{0.00988067050261},{0.00891299549302},{0.0210105783925},{0.00988067050261},{0.00562487715634},{0.0050739982403},{0.011960921317},{0.00891299549302},{0.0050739982403},{0.00457707029451},{0.0107895145135},{0.0210105783925},{0.011960921317},{0.0107895145135},{0.0254340912301},{0.0167842764586},{0.00955496828002},{0.00861919129812},{0.020317994617},{0.00955496828002},{0.00543946109665},{0.00490674111904},{0.0115666465908},{0.00861919129812},{0.00490674111904},{0.00442619369484},{0.0104338535432},{0.020317994617},{0.0115666465908},{0.0104338535432},{0.0245956926577},{0.0157755233003},{0.0335455846742},{0.0157755233003},{0.00898070435783},{0.00810116858408},{0.0190968611776},{0.00810116858408},{0.00730777117393},{0.0172265877666},{0.0335455846742},{0.0190968611776},{0.0172265877666},{0.0406081853165}};

F3x4_data = {{0.035675220653},{0.0183128322366},{0.0182344100227},{0.0382032745233},{0.0244753227002},{0.0125636918383},{0.0125098895363},{0.0262097179792},{0.0243136702555},{0.0124807122786},{0.0124272653253},{0.0260366103541},{0.0443593922569},{0.0227705980131},{0.0226730860234},{0.0475028327522},{0.0169076088239},{0.00867902701785},{0.00864186026479},{0.0181057330441},{0.0115996250193},{0.00595432860953},{0.00592883000691},{0.0124216094776},{0.0115230128428},{0.00591500198703},{0.00588967179532},{0.0123395683309},{0.0210233108084},{0.0107917023874},{0.0107454883894},{0.0225130861001},{0.0175659909532},{0.00901698824865},{0.00897837422257},{0.0188107701193},{0.0120513143089},{0.00618619010983},{0.00615969859186},{0.0129053068343},{0.0119717188549},{0.00614533210901},{0.00611901555987},{0.0128200710061},{0.0218419583344},{0.0112119311773},{0.0111639176066},{0.023389745462},{0.0123371174133},{0.0257370502431},{0.0164887072616},{0.00846399613946},{0.00842775022689},{0.0176571468521},{0.00840809388701},{0.00837208736824},{0.0175405264916},{0.0298843468659},{0.0153402563639},{0.015274563802},{0.0320020419319}};

CF3x4_true = {{0.0416568162879},{0.0196906446876},{0.020059581998},{0.0418708750064},{0.0196906731167},{0.00930752953659},{0.00948192174002},{0.0197918560833},{0.0200595760631},{0.00948190524483},{0.00965956467012},{0.0201626546834},{0.0418709213318},{0.0197918494056},{0.0201626829564},{0.0420860802556},{0.0203517209587},{0.00961999840335},{0.00980024523598},{0.0204563008017},{0.00962001229257},{0.0045472568675},{0.00463245736475},{0.0096694459192},{0.00980024233646},{0.00463244930591},{0.00471924602663},{0.00985060209753},{0.0204563234343},{0.00966944265675},{0.00985061591054},{0.0205614407901},{0.0196808279502},{0.00930287585221},{0.00947718086105},{0.0197819603262},{0.00930288928358},{0.00439735687396},{0.00447974874299},{0.00935069333428},{0.00947717805711},{0.00447974094981},{0.00456367642291},{0.0095258777123},{0.0197819822127},{0.00935069017938},{0.00952589106997},{0.0198836343825},{0.0197967820372},{0.0420965691758},{0.0197968106195},{0.00935769937775},{0.00953303159745},{0.0198985389869},{0.00953301501334},{0.00971163206602},{0.0202713362813},{0.0420966157509},{0.0198985322732},{0.0202713647068},{0.0423129344335}};

CF3x4_data = {{0.0412853327333},{0.0179029185228},{0.0196879766474},{0.0373481379277},{0.023744791903},{0.0102966609843},{0.0113233169635},{0.0214803588671},{0.0264876427489},{0.0114860672931},{0.0126313161929},{0.0239616364764},{0.0430353509497},{0.0186617941686},{0.0205225180086},{0.0389312648494},{0.0195664255595},{0.00848475958363},{0.00933075511284},{0.0177004643578},{0.0112534082309},{0.00487991345405},{0.0053664782087},{0.0101802217625},{0.0125533320378},{0.00544361073975},{0.00598638043206},{0.0113561777357},{0.0203958145674},{0.00884441476503},{0.00972627067101},{0.0184507583003},{0.0203283510097},{0.00881515995474},{0.00969409893203},{0.0183897284359},{0.011691621031},{0.00506993948842},{0.00557545129441},{0.0105766441967},{0.0130421644582},{0.00565558740106},{0.00621949279039},{0.0117983924268},{0.0211900368003},{0.00918882027137},{0.0101050160447},{0.0191692391635},{0.0156818684709},{0.0327146987722},{0.0207989944779},{0.00901924918222},{0.00991853741894},{0.0188154887726},{0.0100610968156},{0.0110642652425},{0.0209889371442},{0.0376963517058},{0.0163465974116},{0.017976478399},{0.0341014217318}};

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


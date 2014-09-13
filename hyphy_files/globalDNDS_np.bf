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

F61_true = {{0.0159744555101},{0.0162616273648},{0.0162584125995},{0.015950404186},{0.0162616273648},{0.0165554372266},{0.0165599318572},{0.0162499968371},{0.0162584125995},{0.0165599318572},{0.0165258239283},{0.0162427541962},{0.015950404186},{0.0162499968371},{0.0162427541962},{0.0159498055864},{0.0159810637454},{0.0165160769394},{0.0162769960766},{0.0161666474214},{0.0165160769394},{0.0168416529064},{0.016823753332},{0.0165250436193},{0.0162769960766},{0.016823753332},{0.0166148969851},{0.0164766049448},{0.0161666474214},{0.0165250436193},{0.0164766049448},{0.0162143431321},{0.016425557835},{0.0165924095857},{0.0167029140915},{0.0162690420722},{0.0165924095857},{0.0168721422871},{0.0168940590607},{0.016559997581},{0.0167029140915},{0.0168940590607},{0.0168842653373},{0.0165622079158},{0.0162690420722},{0.016559997581},{0.0165622079158},{0.0162529412102},{0.0163380426741},{0.0158455302574},{0.0163380426741},{0.0165636807244},{0.0166423313189},{0.0162419356376},{0.0166423313189},{0.0161273090386},{0.0161881807798},{0.0158455302574},{0.0162419356376},{0.0161881807798},{0.0159268238491}};

F61_data = {{0.0173880423911},{0.0164372409044},{0.0176845917926},{0.016179443123},{0.0165702162599},{0.0168751388119},{0.0168612911336},{0.0165827548563},{0.0199979635366},{0.0180027537136},{0.0203125372928},{0.0176880736391},{0.0163758084812},{0.0167143522038},{0.0173526899538},{0.0164437371434},{0.0156511101308},{0.0130772197091},{0.0159224609163},{0.0128498732629},{0.0106227264548},{0.0108344522986},{0.01081450703},{0.0106444749742},{0.0197901564963},{0.0202601728735},{0.0202050806618},{0.0199553363572},{0.0171833970119},{0.0175374850869},{0.0175077803752},{0.0172743615653},{0.0157897439053},{0.013556624406},{0.0160378523109},{0.0133340555461},{0.0196892065367},{0.0200289084326},{0.0200352555623},{0.019673148356},{0.0185637905668},{0.0188040208476},{0.0187655715723},{0.0184725385329},{0.0156897236367},{0.0159838890203},{0.0159669587933},{0.0157103221199},{0.0125903322176},{0.0122313072727},{0.0178324018031},{0.0180708983565},{0.0181472935433},{0.0177427793393},{0.0138471019339},{0.0108928661228},{0.0135039399668},{0.0168857823934},{0.0137964123905},{0.0172136423613},{0.0135464017113}};

F1x4_true = {{0.0141606238625},{0.0157075921214},{0.0150711811496},{0.0144616731666},{0.0157075921214},{0.0174235579342},{0.0167176226545},{0.0160415293634},{0.0150711811496},{0.0167176226545},{0.0160402891462},{0.0153915885442},{0.0144616731666},{0.0160415293634},{0.0153915885442},{0.0147691226607},{0.0157075921214},{0.0174235579342},{0.0167176226545},{0.0160415293634},{0.0174235579342},{0.0193269833301},{0.0185439286042},{0.0177939759357},{0.0167176226545},{0.0185439286042},{0.0177926002317},{0.0170730327491},{0.0160415293634},{0.0177939759357},{0.0170730327491},{0.01638256598},{0.0150711811496},{0.0167176226545},{0.0160402891462},{0.0153915885442},{0.0167176226545},{0.0185439286042},{0.0177926002317},{0.0170730327491},{0.0160402891462},{0.0177926002317},{0.0170717127833},{0.0163812993962},{0.0153915885442},{0.0170730327491},{0.0163812993962},{0.0157188076741},{0.0160415293634},{0.0147691226607},{0.0160415293634},{0.0177939759357},{0.0170730327491},{0.01638256598},{0.0170730327491},{0.0163812993962},{0.0157188076741},{0.0147691226607},{0.01638256598},{0.0157188076741},{0.0150831084104}};

F1x4_data = {{0.0147191558864},{0.0156001275931},{0.01623648803},{0.0144161224981},{0.0156001275931},{0.0165338272657},{0.0172082751815},{0.0152789570342},{0.01623648803},{0.0172082751815},{0.0179102351781},{0.015902216281},{0.0144161224981},{0.0152789570342},{0.015902216281},{0.0141193278666},{0.0156001275931},{0.0165338272657},{0.0172082751815},{0.0152789570342},{0.0165338272657},{0.0175234107811},{0.0182382257897},{0.0161934339893},{0.0172082751815},{0.0182382257897},{0.0189821995336},{0.0168539965819},{0.0152789570342},{0.0161934339893},{0.0168539965819},{0.0149643986346},{0.01623648803},{0.0172082751815},{0.0179102351781},{0.015902216281},{0.0172082751815},{0.0182382257897},{0.0189821995336},{0.0168539965819},{0.0179102351781},{0.0189821995336},{0.019756521456},{0.0175415048451},{0.015902216281},{0.0168539965819},{0.0175415048451},{0.0155748264144},{0.0152789570342},{0.0141193278666},{0.0152789570342},{0.0161934339893},{0.0168539965819},{0.0149643986346},{0.0168539965819},{0.0175415048451},{0.0155748264144},{0.0141193278666},{0.0149643986346},{0.0155748264144},{0.0138286435504}};

F3x4_true = {{0.0139509930507},{0.0162487180251},{0.0151906998939},{0.0159166873802},{0.0162487180251},{0.0189248777129},{0.0176926042671},{0.0185381617061},{0.0151906998939},{0.0176926042671},{0.0165405690067},{0.0173310688651},{0.0159166873802},{0.0185381617061},{0.0173310688651},{0.0181593479574},{0.0141210767902},{0.0164468145128},{0.0153758975378},{0.0161107359114},{0.0164468145128},{0.0191556006413},{0.0179083039155},{0.0187641699805},{0.0153758975378},{0.0179083039155},{0.0167422236},{0.0175423608492},{0.0161107359114},{0.0187641699805},{0.0175423608492},{0.0183807379184},{0.014248432894},{0.0165951461343},{0.0155145707023},{0.0162560364849},{0.0165951461343},{0.0193283624428},{0.0180698165145},{0.0189334015212},{0.0155145707023},{0.0180698165145},{0.0168932194765},{0.0177005730566},{0.0162560364849},{0.0189334015212},{0.0177005730566},{0.0185465113366},{0.0131919478809},{0.0129223800938},{0.0131919478809},{0.0153646583106},{0.0143642047951},{0.0150506927782},{0.0143642047951},{0.0134288947546},{0.014070682797},{0.0129223800938},{0.0150506927782},{0.014070682797},{0.0147431429012}};

F3x4_data = {{0.0142519897201},{0.0153529026809},{0.0151914430195},{0.0150784095804},{0.0178227082998},{0.0191994459308},{0.0189975338817},{0.018856180846},{0.0183714337539},{0.0197905583761},{0.0195824298598},{0.0194367248476},{0.017833448029},{0.0192110152639},{0.0190089815452},{0.0188675433321},{0.0128478822995},{0.0138403332077},{0.0136947805679},{0.0135928831943},{0.0160668133356},{0.0173079146405},{0.0171258949862},{0.0169984680654},{0.0165614783042},{0.0178407906299},{0.0176531669554},{0.0175218168152},{0.0160764949856},{0.0173183441618},{0.0171362148248},{0.0170087111183},{0.0141818756196},{0.015277372528},{0.0151167071837},{0.015004229824},{0.0177350276892},{0.0191049923205},{0.0189040735981},{0.0187634159629},{0.018281053633},{0.0196931967285},{0.0194860921215},{0.0193411039197},{0.0177457145832},{0.0191165047371},{0.0189154649437},{0.0187747225501},{0.010861819815},{0.0106676223619},{0.0126091495655},{0.0135831592619},{0.0134403112063},{0.0133403072373},{0.0140013575013},{0.0138541114405},{0.0137510285498},{0.0126167476728},{0.0135913442947},{0.0134484101607},{0.0133483459306}};

CF3x4_true = {{0.0159551157803},{0.0162859953563},{0.0162197212968},{0.0159532203971},{0.0162859187359},{0.016623658553},{0.0165560104104},{0.0162839840551},{0.0162197001014},{0.0165560666665},{0.0164886935816},{0.0162177732869},{0.0159531468116},{0.0162839855547},{0.0162177196739},{0.0159512516622},{0.0161496307232},{0.0164845441792},{0.0164174621471},{0.0161477122327},{0.0164844666247},{0.0168263239576},{0.0167578510905},{0.0164825083574},{0.0164174406933},{0.0167579080325},{0.0166897135764},{0.0164154903884},{0.01614763775},{0.0164825098754},{0.0164154361217},{0.0161457194962},{0.0162952628822},{0.0166331964796},{0.0165655095234},{0.0162933270913},{0.0166331182258},{0.0169780583178},{0.0169089679845},{0.0166311422995},{0.0165654878762},{0.0169090254399},{0.0168402160283},{0.016563519984},{0.016293251937},{0.0166311438311},{0.016563465228},{0.016291316385},{0.0162453874856},{0.0159134422751},{0.0162453110563},{0.0165822087453},{0.0165147292781},{0.0162433811995},{0.016514785394},{0.0164475802987},{0.0161773355228},{0.015913368873},{0.0162433826954},{0.0161772820434},{0.0159114784491}};

CF3x4_data = {{0.0161917526885},{0.0154265311478},{0.0161074582697},{0.0151507062293},{0.0178456857825},{0.0170022994345},{0.0177527809723},{0.0166982999279},{0.0195585959592},{0.0186342575493},{0.0194567737223},{0.0183010787859},{0.017856455163},{0.0170125598545},{0.0177634942873},{0.0167083768926},{0.014596496685},{0.013906666875},{0.0145205071842},{0.0136580169859},{0.0160874797422},{0.0153271861365},{0.0160037281694},{0.0150531374974},{0.0176316292976},{0.0167983592499},{0.0175398388675},{0.0164980061746},{0.0160971880937},{0.0153364356718},{0.0160133859791},{0.0150622216518},{0.0161120534484},{0.0153505984904},{0.0160281739447},{0.0150761312406},{0.0177578455331},{0.0169186105114},{0.0176653980201},{0.0166161073549},{0.0194623244026},{0.0185425357824},{0.0193610033564},{0.0182109969956},{0.0177685619044},{0.0169288204275},{0.0176760586018},{0.0166261347187},{0.0134225292032},{0.01318253565},{0.0155274206672},{0.0147935954295},{0.0154465847673},{0.0145290873418},{0.0162135520772},{0.0169292183049},{0.0159236552989},{0.0155367910385},{0.0148025229574},{0.0154559063563},{0.0145378552463}};

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


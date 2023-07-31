
# Ecological Niche Modelling with R

This workflow strongly simplifies the process of estimating the spatial probability distribution of species presence/absence given a set of input environmental parameters. It internally embeds 4 models: **AquaMaps (Envelope Model), MaxEnt, Support Vector Machines, Artificial Neural Networks**. It also produces an ensemble model based on the combination of the 4 models and a Biodiversity Index (a species-richness index) by combining the ensemble models of different species. 

It is particularly meant for quick species distribution estimation to be integrated in ecosystem models.

_Input_: environmental parameters in ASC (ESRI-GRID) format, one or more CSV files with species observation indications;

_Output_: Species distributions per individual model (ASC format), one model-ensemble distribution per species (ASC format), one spatial Biodiversity Index dataset (ASC format) for all species.

_Citation_: 

#### Coro G., Sana L., Bove P. An Open Science Automatic Workflow for Multi-model Species Ecological Niche Estimation. https://github.com/cybprojects65/EcologicalNicheModellingWithR

_Auxiliary data and information_:

Sample global-scale environmental parameters are given for 2019. Presence data are reported in CSV format in the Presence Records folders. All species data in this folder will automatically be processed and merged in the Biodiversity Index. 

A Web Processing Service version is available after free subscription to the D4Science e-Infrastructure (BiodiversityLab Virtual Research Environment):

https://services.d4science.org/group/biodiversitylab/data-miner?OperatorId=org.gcube.dataanalysis.wps.statisticalmanager.synchserver.mappedclasses.transducerers.ECOLOGICAL_NICHE_MODELLER

Example input datasets are available at 

https://data.d4science.net/cGHy (environmental parameters)

and 

https://data.d4science.net/MQtR (species occurrence data)

as indicated in the Web interface.

_Instructions_:

The configuration file (_Workflow_Parameters.txt_) commands the entire workflow. In the following, we report the main parameters:

    #environmental data folder containing ASC files:
    env_data_folder="./Environmental Parameters/2019"
    
    #species occurre data folder containing CSV files (see the example in the Presence Record folder for the format):
    presence_data_folder="./Presence Records/"
    
    #output ENM folder:
    output_folder="./ENM"
    
    #environmental data projection folder to project the ENMs on other areas or climatic scenarios
    projection_environmental_layers="./Environmental Parameters 2050 RCP85/"
    
    #output ENM projection folder:
    output_folder_projection="./ENM_p"
    
    #flag to activate/deactivate ecological niche modelling:
    EcologicalNicheModeler =T
    
    #flag to activate/deactivate ensemble model production:
    EnsambleModeler =T
    
    #flag to activate/deactivate biodiversity index production:
    AgreementIndex =T
    
    #number of background point locations to sample:
    abs_sample =10000 
    
    #flag to activate/deactivate native distribution adjustment:
    nativedistribution =T
    
    #flag to activate/deactivate model projection on other areas or climatic scenarios:
    model_projection =F
    
    #flag to activate/deactivate Support Vector Machines:
    SVM_Active =T
    
    #flag to activate/deactivate Artificial Neural Networks:
    ANN_Active =T
    
    #flag to activate/deactivate AquaMaps (Envelope model):
    AQUAMAPS_Active =T
    
    #flag to activate/deactivate MaxEnt:
    MAXENT_Active =T
    
    ###Artificial Neural Network Parameters###
    
    #Number of repetitions for the training
    rp =10
    #threshold for minimum decrease in overall error, default 0.01 = 1%
    thld =0.001
    #maximum steps for the training of the neural network
    stp =1e+06
    
    #training algorithm; possible: backprop, rprop+, rprop-, sag, or slr
    alg =rprop-
    
    #neuron activation function; possible: "logistic" or "tanh" 
    act.fct =logistic
    
    #number of folders for cross-validation
    nfold =1
    
    #number of hidden-layer neurons to test
    hiddens =50,70,100,150,200
    	
    ###Support Vector Machine Parameters###
    
    #costs to test
    costlist =0.001,0.01,0.1,1,10,100
    #kernel functions to test
    kernellistx =linear,polynomial,radial,sigmoid
    
    #coefficients to test in polynomial and sigmoid kernels
    coef0list =0,1
    #gammas to test in polynomial, radial, and sigmoid kernels
    gammalist =0.0001,0.001,0.01,0.1,1,10 
    #degrees to test in polynomial kernels
    degreelist =3,4  
    #number of folders in cross validation
    cross95 =20
    
    ###MaxEnt Parameters###
    #prior species prevalence
    maxent_prevalence =0.5

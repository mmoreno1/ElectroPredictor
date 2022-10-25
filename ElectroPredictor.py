#import libraries
from openbabel import pybel
import numpy as np
import pandas as pd
import subprocess as sp 
import os
import shutil
import weka.core.jvm as jvm 
from sys import exit

#Creating the important files for calculation
def createFiles(molFile):
    '''
    Create 3D structures
    using UFF level of theory 
    
    createFiles(string):
        return: None
    '''
    if molFile.endswith('.sdf'):
        file = pybel.Outputfile('sdf','ToMoCoMD/chemical_datasets/molecules3d.sdf',overwrite = True)
        molecules = pybel.readfile('sdf', molFile)
        for molecule in molecules:
            molecule.make3D(forcefield = 'uff',steps=50)
            molecule.localopt(forcefield = 'uff', steps = 2000)
            file.write(molecule)
        file.close()
    else:
        print("This is not a valid .sdf file")
        print("Please use a valid .sdf file")
    
    return

def getDesc():
    '''
    Get the topographic descriptors for GM
    
    getDesc():
        return pandas DataFrame with topographic descriptors

    '''
    desc = pd.read_csv('ToMoCoMD/Calculation/headings/molecules3d.sdf/molecules3d.sdf_user_specified_headings.csv')
    desc = desc.dropna()
    desc = desc.drop('molecules', axis = 1)
    return desc

def getTCMopac():

    '''
    Get the titles and charges for the molecules
    that are part of GM based on topographic
    descriptors and write them to a .mop file
    for running electronic descriptor calculation

    getTCMopac():
        return list with titles, list with charges
    '''
    executeTomocomd() #execute descriptor calculations
    #Set empty arra
    titles = []
    charges = []
    #Get the molecules
    desc = pd.read_csv('ToMoCoMD/Calculation/headings/molecules3d.sdf/molecules3d.sdf_user_specified_headings.csv')
    desc = desc.dropna()
    inMolecules = desc['molecules'].tolist()
    prefix = 'molecules3d.sdf_'

    for i in range(len(inMolecules)):
        inMolecules[i] = inMolecules[i].removeprefix(prefix)
    #Transform to MOP
    file = open('StructuralOutputFiles/molecules3d.mop','w')
    conv = pybel.ob.OBConversion()
    conv.SetOutFormat("mopcrt")
    conv.AddOption("k",conv.OUTOPTIONS,"PM3")
    molecules  = pybel.readfile("sdf",'ToMoCoMD/chemical_datasets/molecules3d.sdf')
    i = 0
    for molecule in molecules:
        comparator = "_".join(molecule.title.split())
        comparator = comparator.replace('-','_')
        if comparator in inMolecules:
            titles.append(comparator)
            charges.append(molecule.charge)
            out = conv.WriteString(molecule.OBMol)
            file.write(out)
            file.write("\n")
        i = i+1
            
    file.close()
    molecules.close()
    
    return titles, charges
    

def getHLW():
    '''
    
    Return the HOMO,LUMO and w from the molecules after PM3 geometry
    optimization
    
    getHLW():
        return list with homo, list with lumo, list with w
    '''
    executeMopac()
    file = open('StructuralOutputFiles/molecules3d.arc','r')
    newlist = [] #Create a list for appending the array

    for line in file:
        newlist.append(line)
    
    #Set an empty array to try to parse orbitals
    orbitalsWeight = []

    for i in range(len(newlist)):
        if "HOMO" in newlist[i].strip():
            for j in range(i,len(newlist)):
                orbitalsWeight.append(newlist[j].strip())
                if  "MOLECULAR" in newlist[j].strip():
                    break

    #Orbitals and molecular weight were grabbed
        
    #Grab only the orbitals

    orbitalsOnly = [] #set an empty array for parsing orbitals

    for i in range(len(orbitalsWeight)):
        if i%2 == 0:
            orbitalsOnly.append(orbitalsWeight[i])
    #Set empty array for parsing homo and lumo     
    homo = []
    lumo = []
    for i in range(len(orbitalsOnly)):
        line = orbitalsOnly[i]
        lineList = line.split()
        if float(lineList[5]) is not np.nan and float(lineList[6])  is not np.nan:
            homo.append(float(lineList[5]))
            lumo.append(float(lineList[6]))
        else:
            homo.append(np.nan)
            lumo.append(np.nan)

    #HOMO and LUMO in ev
    homo = np.array(homo)
    lumo = np.array(lumo)

    #Transfrm to hartree
    homo = homo/27.211396641308

    lumo = lumo/27.211396641308
    
    #Create an array
    w = (np.power((lumo + homo),2))/(8*(lumo-homo))
    
    homo = homo.tolist()
    lumo = lumo.tolist()
    w = w.tolist()
    
    file.close()
    return homo,lumo,w

def getTCHLW():
    '''

    Concat into a dataframe the titles, charges MOPAC
    descriptors for the molecules

    getTCHLW():
        pandas dataframe with titles, charges, homo, lumo and w
    '''
    titles,charges = getTCMopac()
    homo,lumo,w = getHLW()
    
    labels = ['Titles','Charges','HOMO','LUMO','w']
    tchlw = pd.DataFrame(list(zip(titles,charges,homo,lumo,w)),columns=labels)
    
    return tchlw

def getTCHLWDesc():
    '''

    Concat into a unique dataframe  the dataMop and the
    descriptors
    
    joinTCHLWDesc():
        return dataframe with titles, charges, MOPAC desc, topographicdesc

    '''
    tchlw = getTCHLW()
    desc = getDesc()
    dataFull = pd.concat([tchlw,desc], axis = 1, join = 'inner')
    
    return dataFull


def divideCharges():

    '''
    Create individual dataframes for each model
    part of GM and return a models list depending 
    if there are charged or neutral molecules in the
    dataset

    divideCharges():
        return valid models list
    '''
    dataFull = getTCHLWDesc()
    dataFull = dataFull.dropna() #remove nan molecules
    charged = dataFull.loc[dataFull['Charges'] > 0]
    noncharged = dataFull.loc[dataFull['Charges'] == 0]
    models = []
    if len(charged.index) == 0:
        models = ['M8_LR','M8_RF']
    elif len(noncharged.index) == 0:
        models = ['M7']
    else:
        models = ['M7','M8_LR','M8_RF']
    
    for model in models:
        m = pd.read_csv(f'Models/{model}_training.csv')
        m = m.drop('E', axis = 1)
        head  = m.columns.to_list()
        if model == 'M7':
            m7input = charged.loc[:,head]
            m7input.to_csv(f'Models/{model}_Input.csv', index = False)
        elif model == 'M8_LR' or model == 'M8_RF':
            m8input = noncharged.loc[:,head]
            m8input.to_csv(f'Models/{model}_Input.csv', index = False)

    return models


#Manipulate files, directories and external commands

def createDirectories():
    '''
    Create the Results and StructuralOutputFiles directories
    if they not exist
    
    createDirectories():
        return None
    '''

    resultsPath = 'Results'
    resultsExist = os.path.exists(resultsPath)
    if not resultsExist:
        os.makedirs(resultsPath)
    structPath = 'StructuralOutputFiles'
    structExist = os.path.exists(structPath)
    if not structExist:
        os.makedirs(structPath)
    return

def getSDFfile():
    '''
    Extract the sdf file from the root path
    and return an string with the filename

    def getSDFfile():
        return sdf file as an string
    '''
    rootPath = os.path.dirname(os.path.abspath(__file__))
    sdfFiles = []
    for file in os.listdir(rootPath):
        if file.endswith(".sdf"):
            sdfFiles.append(file)
    if len(sdfFiles) > 1:
        print('\nThere are more .sdf files in the current path\n')
        print('\n Try again using only one .sdf file')
        exit()
    elif len(sdfFiles) == 0:
        print('\nThere are no .sdf files in the current path\n')
        print('\n Try again using only one .sdf file')
        exit()
    else:
        filename = sdfFiles[0]
    return filename
        

def getTomocomdpath():
    '''
    Get the ToMoCoMD path as an string

    getTomocomdpath():
        return path as an string
    '''
    tomocomdPath = os.path.abspath('ToMoCoMD')
    #Check if downloaded
    file_exists = os.path.exists(tomocomdPath + '/ToMoCoMD-CARDD_CLI.jar')
    if not file_exists:
        print('ToMoCoMD-CARDD_CLI.jar is not in ToMoCoMD folder.\n')
        print('Place a valid copy of ToMoCoMD-CARDD_CLI.jar.\n')
        print('Download at: http://tomocomd.com/software/qubils-midas ')
        exit()
    return tomocomdPath

def getAmbitpath():
    '''
    Get the ambit package path as an string

    getAmbitpath():
        return path as an string

    '''
    ambitPath = os.path.abspath('Ambit')
    file_exists = os.path.exists(ambitPath + '/example-ambit-appdomain-jar-with-dependencies.jar')
    if not file_exists:
        print('example-ambit-appdomain-jar-with-dependencies.jar is not in Ambit folder.\n')
        print('Place a valid copy of example-ambit-appdomain-jar-with-dependencies.jar.\n')
        print('Download at: https://github.com/ideaconsult/appdomain')
        exit()
    return ambitPath

def getStructuralpath():
    '''
    Get the StructuralOutputFiles folder path as an string

    getStructuralpath():
        return path as an string

    '''
    outPath = os.path.abspath('StructuralOutputFiles')
    return outPath


def executeMopac():
    '''
    Execute MOPAC as a shell command
    and print an status for the process

    executeMopac():
        return None
    '''
    outPath = getStructuralpath()
    print('\n\nLoading MOPAC calculations ...')
    process = sp.getstatusoutput('mopac ' + outPath +'/molecules3d.mop')
    print(process[1])
    return

def executeTomocomd():
    '''
    Execute ToMoCoMD using java as a shell command
    and print an status for the process
    
    executeTomocomd():
        return None
    '''
    tomocomdPath = getTomocomdpath()
    print('\n\nLoading ToMoCoMD descriptor calculations ...')
    process = sp.getstatusoutput(f'cd {tomocomdPath};java -jar ToMoCoMD-CARDD_CLI.jar')
    print(process[1])
    return

def deleteFiles():
    '''
    Delete all residual files for a calculation

    deleteFiles():
        return None
    '''

    ambitPath = getAmbitpath()
    #Delete the applicability domain 
    filelist = [ f for f in os.listdir(ambitPath) if f.endswith(".csv") ]
    for f in filelist:
        os.remove(os.path.join(ambitPath, f))
        
    #Delete the Input,toWeka,predictions
    modelsPath = os.path.abspath('Models')
    filelist = [ f for f in os.listdir(modelsPath) if f.endswith("Input.csv") ]
    for f in filelist:
        os.remove(os.path.join(modelsPath, f))
    filelist = [ f for f in os.listdir(modelsPath) if f.endswith("toWeka.csv") ]
    for f in filelist:
        os.remove(os.path.join(modelsPath, f))
    filelist = [ f for f in os.listdir(modelsPath) if f.endswith("predictions.txt") ]
    for f in filelist:
        os.remove(os.path.join(modelsPath, f))
        
    #Delete the molecule files
    structuralPath ='StructuralOutputFiles'
    filelist = [ f for f in os.listdir(structuralPath) if f.startswith("molecules3d") ]
    for f in filelist:
        os.remove(os.path.join(structuralPath, f))
        
    #Delete the molecule 3d and the dir with the calculations
    descriptorsPath = 'ToMoCoMD/chemical_datasets'
    filelist = [ f for f in os.listdir(descriptorsPath) if f.startswith("molecules3d") ]
    for f in filelist:
        os.remove(os.path.join(descriptorsPath, f))
    
    #Delete the headings folder
    headingsPath ='ToMoCoMD/Calculation/headings'
    shutil.rmtree(headingsPath, ignore_errors=True)
    
    return


#Start weka modelling
def executeAmbit():

    '''
    Execute ambit depending on the models that are going
    to be used.

    executeAmbit():
        return list with the valid models
    '''
    analysis = ['RANGE','EUCLIDEAN','DENSITY','CITYBLOCK']
    models = divideCharges()
    ambitPath = getAmbitpath()
    
    print('\nLoading Ambit Discovery calculations ...')
    for model in models:
        training = pd.read_csv(f'Models/{model}_training.csv')
        test = pd.read_csv(f'Models/{model}_Input.csv')
        training = training.drop('E', axis = 1)
        head = training.columns.to_list()
        ambitHead  = []
        ambitString = ''
        ambitHead.append('Titles')
        for j in range(len(head)):
            if head[j] != 'Titles':
                ambitHead.append(f'Desc{j}')
                ambitString = ambitString + f'Desc{j}'
                if j != len(head)-1:
                    ambitString = ambitString + ','
        training.columns = ambitHead
        test.columns = ambitHead
        training.to_csv(f'Ambit/{model}_training_AmbitInput.csv',index  = False)
        test.to_csv(f'Ambit/{model}_test_AmbitInput.csv', index = False)
        for technique in analysis:
            command = f'cd {ambitPath};java -jar example-ambit-appdomain-jar-with-dependencies.jar  -m _mode{technique}  -t {model}_training_AmbitInput.csv -s {model}_test_AmbitInput.csv -f {ambitString} -o {model}_{technique}.csv'
            sp.getstatusoutput(command)
    print('Application domain calculations finished')
    
    return models

def makeAmbitConsensus():

    '''
    Create ambit consensus files created when ambit is executed
    the domain is evaluated using 4 techniques and it leaves the
    domain of the model if a molecules is false in more than 2 techniques
    Separated csv are created for each valid model that is going to WEKA

    makeAmbitConsensus():
        return list with the valid models
    '''
    np.random.seed()
    analysis = ['RANGE','EUCLIDEAN','DENSITY','CITYBLOCK']
    models = executeAmbit()
    for model in models:
        test = pd.read_csv(f'Models/{model}_Input.csv')
        domainList = []
        domainList.append(test)
        for technique in analysis:
            results = pd.read_csv(f'Ambit/{model}_{technique}.csv')
            domain = results.iloc[:,-2]
            domainList.append(domain)
        consensus = pd.concat(domainList, axis = 1)
        consensus['Domain Count'] = consensus.iloc[:,-4] + consensus.iloc[:,-3] 
        + consensus.iloc[:,-2] + consensus.iloc[:,-1]
        inDomain = consensus.loc[consensus['Domain Count'] < 2]
        inDomain = inDomain.iloc[:,:-5]
        inDomain['E'] = np.random.randint(-5, 5, inDomain.shape[0])
        inDomain.to_csv(f'Ambit/{model}_inDomain.csv', index = False)
    
    exception = 0 #handle cases in which no molecules are in domain for M7 and M8
    if 'M7' in models:
        m7 = pd.read_csv('Ambit/M7_inDomain.csv')
        if len(m7.index) == 0:
            print('There are 0 molecules in domain for M7\n')
            exception = exception + 1
            if 'M8_LR' not in models or 'M8_RF' not in models:
                print('There are 0 molecules in domain for M7 and 0 noncharged molecules.\n')
                print('Thanks for using ElectroPredictor V1.0')
                deleteFiles()
                exit()
        m7.to_csv('Models/M7_toWeka.csv', index = False)
            
    if 'M8_LR' in models or 'M8_RF' in models:
        #Consensus by crossing both M8 models domain
        m8LR = pd.read_csv('Ambit/M8_LR_inDomain.csv')
        m8RF = pd.read_csv('Ambit/M8_RF_inDomain.csv')
        m8LRinDomain = m8LR['Titles'].values.tolist()
        m8RFinDomain = m8RF['Titles'].values.tolist()
        inGlobalDomain = [value for value in m8LRinDomain if value in m8RFinDomain]
        if len(inGlobalDomain) == 0:
            print('There are 0 molecules in domain for M8\n')
            exception = exception + 1
            if exception == 2:
                print('There are 0 molecules in domain for M7 and M8\n')
                print('Thanks for using ElectroPredictor V1.0')
                deleteFiles()
                exit()
            if 'M7' not in models:
                print('There are 0 molecules in domain for M8 and 0 charged molecules.\n')
                print('Thanks for using ElectroPredictor V1.0')
                deleteFiles()
                exit()

        m8LR = m8LR.set_index('Titles')
        m8RF = m8RF.set_index('Titles')
        m8LR = m8LR.filter(items = inGlobalDomain, axis=0)
        m8RF = m8RF.filter(items = inGlobalDomain, axis=0)
        m8LR.to_csv('Models/M8_LR_toWeka.csv')
        m8RF.to_csv('Models/M8_RF_toWeka.csv')
    
    return models

def makePredictions():
    '''
    Make predictions for the valid models and store them in txts

    makePredictions()
        return list with valid models
    '''
    jvm.start(system_cp=True, packages=True)
    from weka.core.converters import load_any_file
    from weka.classifiers import Classifier,Evaluation,PredictionOutput
    
    models = makeAmbitConsensus()
    for model in models:
        test = load_any_file(f'Models/{model}_toWeka.csv',class_index = 'last')
        test.delete_first_attribute()
        tupleModel = Classifier.deserialize(f'Models/{model}.model')
        prediction_test = PredictionOutput(classname='weka.classifiers.evaluation.output.prediction.PlainText')
        classifier = tupleModel[0]
        evl = Evaluation(test)
        evl.test_model(classifier, test, prediction_test)
        predictions = prediction_test.buffer_content()
        with open(f'Models/{model}_predictions.txt', 'w') as file:
            file.write(predictions)
    return models

def makeEnssemble():
    '''
    Create the enssemble for M8 if there are neutral molecules 
    in the molecules set provided

    makeEnssemble()
        return list with valid models
    '''

    models = makePredictions()
    
    if 'M8_LR' in models or 'M8_RF' in models:
        modelsEnssemble = ['M8_LR','M8_RF']
        enssemble = []
        for model in modelsEnssemble:
            file = open(f'Models/{model}_predictions.txt','r')
            predicted = []
            for line in file:
                splited = line.split()
                predicted.append(float(splited[2]))
             
            enssemble.append(predicted)  
            file.close()    

        enssembleDf = pd.DataFrame(list(zip(enssemble[0],enssemble[1])),columns=modelsEnssemble)
        enssembleDf['E'] = (enssembleDf['M8_LR'] + enssembleDf['M8_RF'])/2
        enssembleDf.to_csv('Models/M8_ENSS_toWeka.csv', index = False)
    
        
    return models

def makeEnssemblePredictions():
    '''
    Make the predictions from neutral molecules using the vote
    classifier from weka

    makeEnssemblePredictions():
        return list with valid models M7 or M8_ENSS
    '''

    from weka.core.converters import load_any_file
    from weka.classifiers import Classifier,Evaluation,PredictionOutput
    
    models = makeEnssemble()
    if 'M8_LR' in models or 'M8_RF' in models:
        model = 'M8_ENSS'
        test = load_any_file(f'Models/{model}_toWeka.csv',class_index = 'last')
        tupleModel = Classifier.deserialize(f'Models/{model}.model')
        prediction_test = PredictionOutput(classname='weka.classifiers.evaluation.output.prediction.PlainText')
        classifier = tupleModel[0]
        evl = Evaluation(test)
        evl.test_model(classifier, test, prediction_test)
        predictions = prediction_test.buffer_content()
        with open(f'Models/{model}_predictions.txt', 'w') as file:
            file.write(predictions)
    
    finalmodels = []
    if 'M7' in models:
        finalmodels.append('M7')
        if 'M8_LR' in models or 'M8_RF' in models:
            finalmodels.append('M8_ENSS')
    else:
        finalmodels.append('M8_ENSS')
            
        
    return finalmodels

def predictElectrophilicity(molFile):
    '''
    Make the electrophilicity predictions that are part of the AD,
    create the csv files with the results, with the titiles, type, and 
    electrophilicty
    
    predictElectrophilicy(molFile):
        return None
    '''

    models = makeEnssemblePredictions()
    if 'M7' in models:
        titlesm7 = pd.read_csv('Models/M7_toWeka.csv')
        titlesm7 = titlesm7['Titles']
        if 'M8_ENSS' in models:
            titlesm8 = pd.read_csv('Models/M8_LR_toWeka.csv')
            titlesm8 = titlesm8['Titles']
            titles = pd.concat([titlesm7,titlesm8]).tolist()
        else:
            titles = titlesm7.tolist()
    else:
        titlesm8 = pd.read_csv('Models/M8_LR_toWeka.csv')
        titlesm8 = titlesm8['Titles']
        titles = titlesm8.tolist()
        
    labels = ['Titles','Type','E_predicted']
    predictions = []
    chargedMolecules = []
    nonchargedMolecules = []
    for model in models:
        file = open(f'Models/{model}_predictions.txt','r')
        for line in file:
            splited = line.split()
            predictions.append(float(splited[2]))
            if model == 'M7':
                chargedMolecules.append('Charged')
            elif model == 'M8_ENSS':
                nonchargedMolecules.append('Neutral')
        file.close()

    chargesColumn = chargedMolecules + nonchargedMolecules
    
    electrophilicity = pd.DataFrame(list(zip(titles,chargesColumn,predictions)),columns=labels)
    electrophilicity = electrophilicity.set_index('Titles')
    outMolecules = electrophilicity.shape[0]
    electrophilicity.to_csv(f'Results/{molFile}.csv', index = True)
    
    
    inMolecules = 0
    molecules = pybel.readfile('sdf',molFile)
    for molecule in molecules:
        inMolecules = inMolecules + 1

    print('\n\n Succesfull electrophilicity calculation')
    print(f'\n InitialMolecules = {inMolecules}  PredictedMolecules = {outMolecules}')
    print('\n Results .csv in Calculations folder')

    return
    
    
    
    

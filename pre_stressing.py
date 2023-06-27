from abaqus import *
from abaqusConstants import *
import displayGroupOdbToolset as dgo
import time
import os
from shutil import copyfile, rmtree
import copy
import numpy as np
from scipy.spatial import KDTree
from functools import partial
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter


def job_submit(
    jobName,
    SDVINI=True, # if SDVINI subroutine is used
    nodalOutputPrecision=SINGLE, # or FULL
    subroutineFile='', # if having subroutine, assign its name, e.g., "subroutine.for".
    numCpus=1, # avoid multiprocessing when using common blocks if Fortran.
    numGPUs=0,
    intruptWithError=False,
    ensureOdbIsClosed=True,
    ):
    # the possible output is the job status, e.g., 'ERROR', 'ABORTED', 'TERMINATED', 'COMPLETED'.
    
    odbName = jobName + '.odb'
    staName = jobName + '.sta'
    lckName = jobName + '.lck'
    
    if os.path.exists(lckName):
        raise Exception('Error: an lck file is detected!')
    
    if ensureOdbIsClosed and os.path.exists(odbName):
        odbName = jobName+'.odb'
        
        if odbName in session.odbs:
            close_odb(odbName)
    
    if os.path.exists(staName):
        try:
            os.remove(staName)
        except OSError:
            print 'Warning: failed to delete the sta file (%s)! Trying once more after 30 Sec.'
            time.sleep(30)
            os.remove(staName)
    
    if SDVINI:
        modelName = mdb.jobs[jobName].model
        keywordBlock = mdb.models[modelName].keywordBlock
        keywordBlock.setValues(edited = 0)
        keywordBlock.synchVersions(storeNodesAndElements=False)
        
        def _give_idx(searchCriteria):
            return [line.startswith(searchCriteria) for line in keywordBlock.sieBlocks].index(True)
        try:
            idx = _give_idx('** ----------------------------------------'\
                            '------------------------\n** \n** STEP:')
            keywordBlock.insert(idx-1, '\n*INITIAL CONDITIONS, TYPE=SOLUTION, USER')
        except:
            raise Exception(
                'Error: "job_submit" function cannot insert the relavant SDVINI keyboard. '
                'Either fix the function, or insert SDVINI = False, if possible.'
                )
    
    job = mdb.jobs[jobName]
    job.setValues(
        nodalOutputPrecision=nodalOutputPrecision,
        userSubroutine=os.path.join(os.getcwd(), subroutineFile),
        numCpus=numCpus,
        numGPUs=numGPUs
        )
    
    try:
        job.submit(consistencyChecking=ON)
        # job.waitForCompletion() should not be placed just after
        # job.submit(), as it might cause rare strange errors.
        # Instead, the code checks constantly the lck file. 
        # job.waitForCompletion() is still neccessary, otherwise
        # the status of the job in CAE may not be updated sometimes.
        print '%s is just submitted!'%jobName
        
        while not os.path.exists(lckName):
            time.sleep(2)
        
        # printing the sta during submission, while waiting for the lck to be deleted.
        pos = 0L
        
        while os.path.exists(lckName):
            time.sleep(5)
            
            if os.path.exists(staName):
                try:
                    with open(staName, "r") as f:
                        f.seek(pos)
                        
                        for line in f.readlines():
                            print line.strip()
                        
                        pos = f.tell()
                except OSError as ex:
                    print 'Warning: an os error was catched (%s)!' % ex
                finally:
                    f.close()
        
        job.waitForCompletion()
        status = str(job.status)
        
    except AbaqusException as ex:
        job.waitForCompletion()
        
        if intruptWithError == True:
            raise ex
        else:
            print 'Warning: an error is avoided during job submitting (i.e., %s)' % ex
            
            while os.path.exists(lckName):
                print 'Waiting for the lck file to be deleted!'
                time.sleep(10)
            
            status = 'ERROR'
    
    time.sleep(5)
    
    return status


def copy_model(NameOfNewModel, nameOfCopiedModel):
    mdb.Model(name = NameOfNewModel, objectToCopy = mdb.models[nameOfCopiedModel])


def copy_job(NameOfNewJob, nameOfCopiedJob):
    mdb.Job(name = NameOfNewJob, objectToCopy = mdb.jobs[nameOfCopiedJob])


def open_odb(odbName, readOnly=True):
    _open_odb_fn = lambda: session.openOdb(name=odbName, readOnly=readOnly)
    try:
        return _open_odb_fn()
    except Exception as ex1: # if it does not work try again after 5Sec.
        try:
            return _open_odb_fn()
        except Exception as ex2:
            raise Exception('Error: open_odb() did not work. \n%s\n%s'%(ex1, ex2))


def close_odb(odbName, saveOdb = False):
    try:
        odbObj = session.odbs[odbName]
        
        if saveOdb == True:
            odbObj.save()
        
        odbObj.close()
    except Exception as ex:
        raise Exception('Error: close_odb() did not work.\n%s'%(ex))


def close_all_odbs():
    for odbObj in session.odbs.values():
        odbObj.close()


def unit_vector(vector):
    norm = np.linalg.norm(vector)
    if (norm==0):
        vector += 1e-6
        norm = np.linalg.norm(vector)
    
    return vector / norm


class nonlipls_tools():
    
    def __init__(
        self,
        jobName,
        modelName,
        materialName='CAR_UMAT',
        sectionName='CAR',
        stepName='EQ',
        txtFolderName='txt',
        subroutineFile='subroutines.for',
        numSdv=200,
        ):
        
        # ensuring the step is not surrpressed.
        modelObj = mdb.models[modelName]
        stepObj = modelObj.steps[stepName]
        
        if stepObj.suppressed:
            stepObj.resume()
        
        # fixing the outputs viriables of the fieldOutputRequest.
        variables = modelObj.steps[stepName].fieldOutputRequestStates[stepName].variables
        variables = set(variables) | {'LE', 'U', 'UR', 'SDV', 'COORD'}
        
        modelObj.fieldOutputRequests[stepName].setValues(variables=tuple(variables))
        
        copy_model(modelName + '-Backup', modelName)
        
        self.odbNameWithoutOptimizaion = jobName +'.odb'
        self.jobNameWithoutOptimizaion = jobName
        self.modelNameWithoutOptimizaion = modelName
        self.modelName = modelName + '-withEQ'
        self.odbName = jobName + '-withEQ.odb'
        self.jobName = jobName + '-withEQ'
        self.txtFolderName = txtFolderName
        self.modelNameInverse = modelName + '-Inverse'
        self.odbNameInverse = jobName + '-Inverse.odb'
        self.jobNameInverse = jobName + '-Inverse'
        self.stepName = stepName
        self.materialName = materialName
        self.numSdv = numSdv
        self.sectionName = sectionName
        
        self.job_submit = partial(
            job_submit, subroutineFile=subroutineFile, intruptWithError=True
            )
    
    def _set_umat_key(self, modelObj, umat_key):
        # ensure correct material and section definitions
        modelObj.fieldOutputRequests[self.stepName].setValues(position=INTEGRATION_POINTS)
        modelObj.sections[self.sectionName].setValues(material=self.materialName, thickness=None)
        materialObj = modelObj.materials[self.materialName]
        materialObj.Depvar(n = self.numSdv)
        userMaterialObj = materialObj.userMaterial
        props = copy.deepcopy(userMaterialObj.mechanicalConstants)
        props[2] = umat_key
        userMaterialObj.setValues(mechanicalConstants=props)
    
    def _focus_on_first_step(self, modelObj):
        
        for key in modelObj.steps.keys():
            stepObj = modelObj.steps[key]
            if key == 'Initial':
                continue
            if key == self.stepName:
                stepObj.resume()
            else:
                stepObj.suppress()
    
    def _resume_all_steps(self, modelObj):
        
        for key in modelObj.steps.keys():
            stepObj = modelObj.steps[key]
            if key != 'Initial':
                stepObj.resume()
    
    def initialize_params(self, *keys):
        # finding the normalized depth and local plan, and initializing the parameters.
        
        if os.path.exists(self.txtFolderName):
            rmtree(self.txtFolderName)
    
        os.makedirs(self.txtFolderName)
        
        with open(os.path.join(self.txtFolderName, 'DATA.txt'), "w") as f:
            f.write('0\n')
        
        modelObj = mdb.models[self.modelNameWithoutOptimizaion]
        self._set_umat_key(modelObj, 0.0)
        self._focus_on_first_step(modelObj)
        
        self.job_submit(jobName=self.jobNameWithoutOptimizaion)
        
        odb = open_odb(self.odbNameWithoutOptimizaion)
        
        assemblyObj = odb.rootAssembly
        frameObj0 = odb.steps[self.stepName].frames[0]
        
        # get the instance name as it is converted by Abaqus to a general instance for all the model in ODB
        instanceOdbName = odb.rootAssembly.instances.keys()[0]
        
        # get all the nodeset region object
        regionNodeSets = odb.rootAssembly.instances[instanceOdbName].nodeSets
        regionElementSets = odb.rootAssembly.instances[instanceOdbName].elementSets
        
        def _get_coords_of_a_nodeset(nodeSetName):
            valueObj = frameObj0.fieldOutputs['COORD'].getSubset(region=regionNodeSets[nodeSetName],
                                                                 position=NODAL).values
            return np.array([i.data for i in valueObj], dtype=np.float32)
        
        def _point_data_of_a_substructure(cartilageKey = 'LAT_CARTILAGE'):
            
            topCoords = _get_coords_of_a_nodeset('TOP_'+cartilageKey)
            bottomCoords = _get_coords_of_a_nodeset('BOTTOM_'+cartilageKey)
            bottomTree = KDTree(bottomCoords)
            topTree = KDTree(topCoords)
            centralPoint = np.mean(np.concatenate((topCoords, bottomCoords)), axis = 0, dtype=np.float32)
            centralPoint += 1e-8 # just to avoid zero devision
            
            ####### extracting integration point data #######
            pointSetCoords = []
            for fieldName in ['SDV91', 'SDV92', 'SDV93']:
                valueObj = frameObj0.fieldOutputs[fieldName].getSubset(region=regionElementSets[cartilageKey],
                                                                       position=INTEGRATION_POINT).values
                if fieldName == 'SDV91': # only once in the loop:
                    coordsFromSdv = np.array([[i.elementLabel, i.integrationPoint, i.data] for i in valueObj],
                                             dtype=np.float32)
                else:
                    coordsFromSdv = np.array([[i.data] for i in valueObj], dtype=np.float32)
                
                pointSetCoords.append(coordsFromSdv)
            
            temp1, temp2, temp3 = np.hsplit(pointSetCoords[0], 3)
            pointSetCoords = np.array([temp1, temp2, temp3, pointSetCoords[1], pointSetCoords[2]],
                                      dtype=np.float32)
            pointSetCoords = np.squeeze(pointSetCoords).T # deleting redundant axis and transposing
            
            def _get_one_unit(point, numPoints = 1):
                '''A helper function:
                   "point" is the array of points,
                   "numPoints" is the number of nodal points averaged'''
                topDistance , topId = topTree.query(point, numPoints)
                bottomDistance , bottomId = bottomTree.query(point, numPoints)
                
                topDistance = np.mean(topDistance)
                bottomDistance = np.mean(bottomDistance)
                
                totalDistance = topDistance + bottomDistance
                depth = topDistance / totalDistance
                
                topPoint = topCoords[topId]
                bottomPoint = bottomCoords[bottomId]
                unit1 = unit_vector(bottomPoint - topPoint)   # unit of the depth (first important unit)
                vec3 = np.cross(unit1, topPoint - centralPoint)   # perpendicular to the surface of central point and depth vector.
                unit3 = unit_vector(vec3)
                vec2 = np.cross(unit3, unit1) # in the surface of central point and depth vector but perpendicular to the depth vector
                unit2 = unit_vector(vec2) # second important unit
                
                return np.hstack([depth, unit1, unit2, unit3])
            
            # Output is element label, integration point lable, depth, units for each point
            return np.array([np.concatenate([pointData[0:2], _get_one_unit(pointData[2:])])
                             for pointData in pointSetCoords], dtype=np.float32)
        
        def _point_sdv_data(point_data):
            '''It takes point_data, i.e., [element, integration point, depth, and units]
               and returns [element, integration point, and relevant sdvs]'''
            
            element = point_data[0]
            integrationPoint = point_data[1]
            depth = point_data[2]
            unit1 = point_data[3:6]
            unit2 = point_data[6:9]
            unit3 = point_data[9:12]
            
            sdv2=1.4 * (depth ** 2) - 1.1 * depth + 0.59 # FIBER CONSTANT
            sdv3=0.1 + 0.2 * depth # SOLID MATERIAL CONSTANT
            alpha1=[None] * 10
            alpha1[0]=0.005
            alpha1[1]=0.01
            alpha1[2]=0.025
            alpha1[3]=0.035
            alpha1[4]=0.042
            alpha1[5]=0.048
            alpha1[6]=0.053
            alpha1[7]=0.058
            alpha1[8]=0.06
            alpha1[9]=0.06  # the deepest points have index 9 but it is still part its upper leyer (with index 8).
            sdv1=alpha1[int(depth*9)] # GAG CONSTANT
            
            # for primary fibril vectors 1 and 2
            pFibrilVec1 = depth*unit1 + (1.0-depth)*unit2
            pFibrilVec2 = depth*unit1 - (1.0-depth)*unit2
            pFibrilVec1 = unit_vector(pFibrilVec1)
            pFibrilVec2 = unit_vector(pFibrilVec2)
            
            sFibrilVec1 = unit1
            sFibrilVec2 = unit2
            sFibrilVec3 = unit3
            sFibrilVec4 = unit_vector( unit1 + unit2 + unit3)
            sFibrilVec5 = unit_vector(-unit1 + unit2 + unit3)
            sFibrilVec6 = unit_vector( unit1 - unit2 + unit3)
            sFibrilVec7 = unit_vector( unit1 + unit2 - unit3)
            
            return np.hstack([
                element, integrationPoint, sdv1, sdv2, sdv3, pFibrilVec1,
                pFibrilVec2, sFibrilVec1, sFibrilVec2, sFibrilVec3, sFibrilVec4,
                sFibrilVec5, sFibrilVec6, sFibrilVec7
                ])
    	
        data = []
        for key in keys:
            subData = _point_data_of_a_substructure(key)
            
            # Depth is normalized within [0,1] (as it is already around 0.1 to 0.9)
            subData[:,2] = (subData[:,2] - np.min(subData[:,2]))/np.ptp(subData[:,2])
            
            data.append(
                np.array([_point_sdv_data(point_data) for point_data in subData], dtype=np.float32)
                )
        
        data = np.concatenate(data, axis = 0)
        
        elementArray = np.unique(data[:,0])
        
        for element in elementArray:
            np.savetxt(os.path.join(self.txtFolderName, '%i.txt'%(element)),
                       data[data[:,0] == element][:,1:],
                       delimiter=',',
                       fmt='%10.7f')
        
        with open(os.path.join(self.txtFolderName, 'DATA.txt'), "w") as f:
            f.write('1\n')
        
        self._resume_all_steps(modelObj)
        self._set_umat_key(modelObj, 1.0)
        
        print 'INITIALIZATION IS COMPLETED! \n'
    
    def run_prestress_optimizer(
            self,
            key,
            sdvList=['SDV%s'%(i) for i in (range(1,4) + range(16,43))],
            zeta=1.0, 
            breakPoint=0, 
            errorLimit=1e-3,
            maxiteration=50,
            eta=4.0,
            ):
        # main pre-stress function
        zeta = float(zeta) # avoiding problem with integer division
        modelWithoutOptimizaionObj = mdb.models[self.modelNameWithoutOptimizaion]
        self._set_umat_key(modelWithoutOptimizaionObj, 1.0)
        self._focus_on_first_step(modelWithoutOptimizaionObj)
        
        # ensure the correct naming of all cartilage sets.
        nodeSet = key + '_NODES'
        elementSet = key + '_ELEMENTS'
        
        # creating helper models, jobs, sets, and BCs
        rootAssemblyObj = modelWithoutOptimizaionObj.rootAssembly
        nodeObj = rootAssemblyObj.sets[nodeSet].nodes
        nodeNameList = ['TEMP-%s'%(i) for i in xrange(1, len(nodeObj)+1)]
        
        for i in xrange(len(nodeNameList)):
            rootAssemblyObj.Set(nodes=(nodeObj[i:i+1],), name=nodeNameList[i])
        
        for name in [self.modelName, self.modelNameInverse]:
            copy_model(name, self.modelNameWithoutOptimizaion)
        
        for name in [self.jobName, self.jobNameInverse]:
            copy_job(name, self.jobNameWithoutOptimizaion)
        
        mdb.jobs[self.jobNameInverse].setValues(model=self.modelNameInverse)
        mdb.jobs[self.jobName].setValues(model=self.modelName)
        steps = mdb.models[self.modelName].steps
        stepsWithoutEQ = modelWithoutOptimizaionObj.steps
        
        BCobj = mdb.models[self.modelNameInverse].boundaryConditions
        BCkeys = BCobj.keys()
        
        for i in BCkeys:
            BCobj[i].suppress()
        
        modelObjTemp = mdb.models[self.modelNameInverse]
        
        for i in xrange(len(nodeNameList)):
            modelObjTemp.VelocityBC(name=nodeNameList[i], # bc name and node names are the same.
                                    createStepName='Initial',
                                    region=modelObjTemp.rootAssembly.sets[nodeNameList[i]],
                                    v1=SET,
                                    v2=SET,
                                    v3=SET,
                                    amplitude=UNSET,
                                    localCsys=None,
                                    distributionType=UNIFORM,
                                    fieldName='')
        
        # some helper functions
        def _integration_points_values(odb, parameters=['SDV3'], frameNum=-1):
            
            instanceOdbName = odb.rootAssembly.instances.keys()[0] # refers to all
            regionElementSets = odb.rootAssembly.instances[instanceOdbName].elementSets
            frameObj = odb.steps[self.stepName].frames[frameNum]
            
            return [np.ravel([[item.data, item.elementLabel, item.integrationPoint]
                    for item in frameObj.fieldOutputs[sdvName].getSubset(region=
                                regionElementSets[elementSet], position=
                                INTEGRATION_POINT).values])
                    for sdvName in parameters]
        
        def _extract_coords_values(odb, frameNum = -1):
            instanceOdbName = odb.rootAssembly.instances.keys()[0] # refers to all
            regionNodeSets = odb.rootAssembly.instances[instanceOdbName].nodeSets
            frameObj = odb.steps[self.stepName].frames[frameNum]
            
            # a list of [[node label], [node coord 1, node coord 2, ...], ...],
            # then flatten the list.
            return [nodeDataElement
                    for nodeName in nodeNameList
                    for nodeData in ([nodeName],
                                     frameObj.fieldOutputs['COORD'].getSubset(region
                                        =regionNodeSets[nodeName], position=NODAL
                                        ).values[0].data.tolist()
                                    )
                    for nodeDataElement in nodeData]
        
        def _edit_node_by_offset(displacementFromInitial, modelName):
            rootAssemblyObj = mdb.models[modelName].rootAssembly
            num = 0
            for i in displacementFromInitial:
                num += 1
                if num % 4 == 1:
                    nodeLabel = i
                elif num % 4 == 2:
                    u1 = -i
                elif num % 4 == 3:
                    u2 = -i
                elif num % 4 == 0:
                    u3 = -i
                    rootAssemblyObj.editNode(nodes=rootAssemblyObj.sets[nodeLabel].nodes[0],
	                                         offset1=u1,
	                                         offset2=u2,
	                                         offset3=u3,
	                                         projectToGeometry=OFF)
        
        def _inverse_run(displacementFromInitial):
            ModelObjTemp = mdb.models[self.modelNameInverse]
            bcStateObj = ModelObjTemp.steps[self.stepName].boundaryConditionStates
            num = 0
            for i in displacementFromInitial:
                num += 1
                if num % 4 == 1:
                    bcLabel = i
                elif num % 4 == 2:
                    # v1 = bcStateObj[bcLabel].v1-i
                    v1 = -i
                elif num % 4 == 3:
                    # v2 = bcStateObj[bcLabel].v2-i
                    v2 = -i
                elif num % 4 == 0:
                    # v3 = bcStateObj[bcLabel].v3-i
                    v3 = -i
                    ModelObjTemp.boundaryConditions[bcLabel].setValuesInStep(stepName=self.stepName,
                                                                             v1=v1,
                                                                             v2=v2,
                                                                             v3=v3)
            
            with open(os.path.join(self.txtFolderName, 'DATA.txt'), "w") as f:
                f.write('-1\n')
            
            self.job_submit(self.jobNameInverse)
            
            with open(os.path.join(self.txtFolderName, 'DATA.txt'), "w") as f:
                f.write('1\n')
        
        def _new_SDV_in_fortran(newSdvData):
            IntegrationPointArray = np.unique(newSdvData[0][2::3]) # [1.0, 2.0, 3.0, 4.0, ...]
            IntegrationCount = IntegrationPointArray[-1].max() # number of all integration points
            _, ElementsIdx = np.unique(newSdvData[0][1::3], return_index=True)  # e.g., [0L, 27L, 54L, ...]
            elementCount = ElementsIdx[1] # all elements have the same number of nodes
            
            for ElementsIdxItem in ElementsIdx:
                elementIdxArray = ElementsIdxItem*3 + 1 + np.arange(0, 3*IntegrationCount, 3, dtype=int)
                ValueIdxArray = elementIdxArray - 1
                IntegrationPointIdxArray = elementIdxArray + 1
                elementItem = newSdvData[0][elementIdxArray[0]]
                sdvDataList = np.concatenate((IntegrationPointArray.reshape(1,-1),
                                              np.take(newSdvData, ValueIdxArray, axis = -1)))
                np.savetxt(os.path.join(self.txtFolderName, '%i.txt'%(elementItem)),
                           sdvDataList.T,
                           delimiter=',',
                           fmt='%10.7f')
        
        if self.job_submit(self.jobNameWithoutOptimizaion) == 'ABORTED':
            raise Exception('ERROR! TOTALLY UNSTABLE MODEL')
        
        odbObjWithoutOptimization = open_odb(self.odbNameWithoutOptimizaion)
        initialNodalCoords = _extract_coords_values(odbObjWithoutOptimization, 0)
        copyfile(self.odbNameWithoutOptimizaion, self.odbName)
        
        def _calculate_r_u(zeta):
            odb = open_odb(self.odbName)
            newNodalCoords = _extract_coords_values(odb, -1)
            close_odb(self.odbName)
            
            displacementFromInitial = [(newNodalCoords[i]-initialNodalCoords[i])*zeta
                                       if i % 4 != 0
                                       else newNodalCoords[i] # just the label
                                       for i in xrange(len(newNodalCoords))]
            
            newError = np.linalg.norm(np.array(
                [displacementFromInitial[i] for i in xrange(len(displacementFromInitial)) if i % 4 != 0]
                )) / zeta
            
            return newError, displacementFromInitial
            
        newError, displacementFromInitial = _calculate_r_u(zeta)
        
        self.optimizerStatus = {'step': [1], 'error': [newError], 'zeta': [zeta]}
        
        failed=False
        iterationNumber = 1
        
        newSdvDataBackup = _integration_points_values(odbObjWithoutOptimization, sdvList, 0)
        close_odb(self.odbNameWithoutOptimizaion)
        previousError = newError
        
        while True:
            
            if iterationNumber == maxiteration:
                failed = True
                break
            else:
                iterationNumber += 1
            
            copy_model(self.modelName + '-Backup', self.modelName)
            _edit_node_by_offset(displacementFromInitial, self.modelName)
            copy_model(self.modelNameInverse + '-Backup', self.modelNameInverse)
            _inverse_run(displacementFromInitial)
            _edit_node_by_offset(displacementFromInitial, self.modelNameInverse)
            odbInverse = open_odb(self.odbNameInverse)
            newSdvData = _integration_points_values(odbInverse, sdvList, -1)
            close_odb(self.odbNameInverse)
            _new_SDV_in_fortran(newSdvData)
            
            
            if self.job_submit(self.jobName) != 'ABORTED':
                newError, displacementFromInitial = _calculate_r_u(zeta)
                if previousError < newError:
                    successfulStep = False
                else:
                    successfulStep = True
                
            else:
                successfulStep = False
            
            print '\n** #STEP: %s | ERROR: %s | ZETA: %s **\n' % (iterationNumber, newError, zeta)
            
            self.optimizerStatus['step'].append(iterationNumber)
            self.optimizerStatus['error'].append(newError)
            self.optimizerStatus['zeta'].append(zeta)
            
            if errorLimit > newError:
                failed = False
                break
            
            if successfulStep == True:
                newSdvDataBackup = copy.deepcopy(newSdvData)
                previousError = newError
            
            else:
                _new_SDV_in_fortran(newSdvDataBackup)
                zeta = zeta/eta
                
                if zeta < 0.0001:
                    failed = True
                    break
                
                copy_model(self.modelName, self.modelName + '-Backup')
                copy_model(self.modelNameInverse, self.modelNameInverse + '-Backup')
        
        # finish_optimization
        modelsObj = mdb.models
        self._resume_all_steps(modelsObj[self.modelName])
        self._resume_all_steps(modelsObj[self.modelNameWithoutOptimizaion])
        
        del modelsObj[self.modelName+'-Backup']
        del modelsObj[self.modelNameInverse]
        del modelsObj[self.modelNameInverse+'-Backup']
        del mdb.jobs[self.jobNameInverse]
        
        for tempSetName in nodeNameList:
            for modelName in [self.modelName, self.modelNameWithoutOptimizaion]:
                del modelsObj[modelName].rootAssembly.sets[tempSetName]
        
        if failed == True:
            print 'PRE_STRESSING HAS NOT BEEN FULLY CONVERGED! \n'
            return 'ABORTED'
        else:
            print 'PRE_STRESSING HAS BEEN COMPLETED! \n'
            return 'COMPLETED'


os.chdir('C:\\Temp\\pre_stress_3d')
openMdb(pathName = 'open_knee.cae')

start_time = time.time()

# run PSA
nonliplsTools = nonlipls_tools('knee', 'knee')
nonliplsTools.initialize_params('LAT_CARTILAGE', 'MED_CARTILAGE', 'FEMUR_CARTILAGE')
nonliplsTools.run_prestress_optimizer('ARTICULAR_CARTILAGE')

elapsed_time = time.time() - start_time
print "Runtime: %i hour(s) and %i minute(s)"%(elapsed_time//3600, (elapsed_time%3600)//60)


# img dir
os.mkdir('img')

# Convergence plot:
SMALL_SIZE = 20//1.4
MEDIUM_SIZE = 24//1.4
BIGGER_SIZE = 28//1.4

plt.rc('font', size=SMALL_SIZE)         
plt.rc('axes', titlesize=MEDIUM_SIZE)    
plt.rc('axes', labelsize=MEDIUM_SIZE)   
plt.rc('xtick', labelsize=SMALL_SIZE)   
plt.rc('ytick', labelsize=SMALL_SIZE)   
plt.rc('legend', fontsize=SMALL_SIZE)   
plt.rc('figure', titlesize=BIGGER_SIZE)

x = nonliplsTools.optimizerStatus['step']
y1 = nonliplsTools.optimizerStatus['error']
y2 = nonliplsTools.optimizerStatus['zeta']

fig, ax1 = plt.subplots()

ax1.plot(x, y1, color='blue', linestyle='-', linewidth=2)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax1.tick_params(axis='y')
ax1.set_xlabel('Step')
ax1.set_ylabel('Error')

ax1.yaxis.set_major_formatter(FuncFormatter(lambda y, _: '{:.2f}'.format(y))) 
ax1.tick_params(axis='y')
fig.tight_layout()
plt.savefig('img/convergence_plot.png', dpi=300)
# plt.show()

# Creating variables to report
viewportObj = session.viewports['Viewport: 1']
session.printOptions.setValues(reduceColors=False)
session.pngOptions.setValues(imageSize=(4000, 2289))


odbObj = open_odb('knee-withEQ.odb')
fieldOutputsObj = odbObj.steps['EQ'].frames[-1].fieldOutputs

s1f1_SDV70 = fieldOutputsObj['SDV70']
s1f1_SDV71 = fieldOutputsObj['SDV71']
s1f1_SDV72 = fieldOutputsObj['SDV72']
s1f1_SDV73 = fieldOutputsObj['SDV73']
s1f1_SDV74 = fieldOutputsObj['SDV74']
s1f1_SDV75 = fieldOutputsObj['SDV75']

scratchOdb = session.ScratchOdb(odb=odbObj)
sessionStep = scratchOdb.Step(name='Session Step', 
    description='Step for Viewer non-persistent fields', 
    domain=TIME, timePeriod=1.0)
sessionFrame = sessionStep.Frame(frameId=0, frameValue=0.0, description='Session Frame')

sessionFrameObj = session.scratchOdbs['knee-withEQ.odb'].steps['Session Step'].frames[0]
sessionFrameObj.FieldOutput(
    name='Fibrillar Mises Stress (MPa)', 
    description='',
    field=sqrt(0.5*((s1f1_SDV70-s1f1_SDV71)*(s1f1_SDV70-s1f1_SDV71)+(
        s1f1_SDV71-s1f1_SDV72)*(s1f1_SDV71-s1f1_SDV72)+(
        s1f1_SDV70-s1f1_SDV72)*(s1f1_SDV70-s1f1_SDV72)+6*(
        s1f1_SDV73*s1f1_SDV73+s1f1_SDV74*s1f1_SDV74+s1f1_SDV75*s1f1_SDV75)))
    )

sessionFrameObj.FieldOutput(
    name='Deformation magnitude (mm)', 
    description='',
    field=fieldOutputsObj['U'].getScalarField(invariant=MAGNITUDE)
    )


def print_png(name):
    session.printToFile(
        fileName='img/'+name+'.png', format=PNG, canvasObjects=(viewportObj, )
        )

leaf_tibia = dgo.LeafFromOdbElementSections(elementSections=(
        'PART-1-1.TIBIA_CARTILAGE_LAT-1__PICKEDSET21.Section-TIBIA_CARTILAGE_LAT-1__PICKEDSET21', 
        'PART-1-1.TIBIA_CARTILAGE_MED-1__PICKEDSET26.Section-TIBIA_CARTILAGE_MED-1__PICKEDSET26', 
        ))
leaf_femur = dgo.LeafFromOdbElementSections(elementSections=(
    'PART-1-1.FEMUR_CARTILAGE-1__PICKEDSET25.Section-FEMUR_CARTILAGE-1__PICKEDSET25', 
    ))
leaf_cartilage = dgo.LeafFromOdbElementMaterials(elementMaterials=('CAR_UMAT', ))

replace_leaf_fn = lambda leaf: viewportObj.odbDisplay.displayGroup.replace(leaf=leaf)

# drawing contours using the following helper functions
# note that these functions are set to my screen; you can adjust them.
def tibia_b(name):
    replace_leaf_fn(leaf_tibia)
    viewportObj.view.setProjection(projection=PERSPECTIVE)
    viewportObj.view.setValues(nearPlane=137.258, 
        farPlane=200.913, width=108.151, height=61.9818, cameraPosition=(
        104.938, -114.507, 65.1538), cameraUpVector=(0.922775, 0.376379, 
        -0.0826106), cameraTarget=(95.8773, 53.5151, 46.2021), 
        viewOffsetX=1.36828, viewOffsetY=-1.61421)
    print_png(name+'_tibia_b')


def femur_b(name):
    replace_leaf_fn(leaf_femur)
    viewportObj.view.setProjection(projection=PERSPECTIVE)
    viewportObj.view.setValues(nearPlane=102.691, 
        farPlane=196.321, width=110.252, height=63.1861, cameraPosition=(
        56.5754, 219.107, 23.7612), cameraUpVector=(0.97178, -0.177, 0.155934), 
        cameraTarget=(80.6936, 52.6705, 43.413), viewOffsetX=1.12, 
        viewOffsetY=-0.107319)
    print_png(name+'_femur_b')


def full_cartilage(name):
    replace_leaf_fn(leaf=leaf_cartilage)
    viewportObj.view.setValues(nearPlane=63.5544, 
        farPlane=186.566, width=72.5957, height=41.605, cameraPosition=(
        178.436, 85.7841, -35.0266), cameraUpVector=(-0.455834, 0.87901, 
        0.139848), cameraTarget=(54.1926, 66.9635, 78.453), 
        viewOffsetX=1.13032, viewOffsetY=-5.57586)
    viewportObj.view.setProjection(projection=PARALLEL)
    viewportObj.view.fitView()
    print_png(name+'_full_cartilage')

viewportObj.setValues(displayedObject=odbObj)
viewportObj.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
viewportObj.viewportAnnotationOptions.setValues(
    legendBox=OFF, triad=OFF, legendNumberFormat=FIXED, title=OFF, state=OFF
    )

odbname = viewportObj.odbDisplay.name
frame1 = session.scratchOdbs[odbname].steps['Session Step'].frames[0]
viewportObj.odbDisplay.setFrame(frame=frame1)

viewportObj.odbDisplay.contourOptions.setValues(
    numIntervals=6, maxAutoCompute=OFF, maxValue=0.12, minAutoCompute=OFF, minValue=0.0
    )
label = 'Fibrillar Mises Stress (MPa)'
viewportObj.odbDisplay.setPrimaryVariable(variableLabel=label, outputPosition=INTEGRATION_POINT)
tibia_b(label)
femur_b(label)

label = 'Deformation magnitude (mm)'
viewportObj.odbDisplay.setPrimaryVariable(variableLabel=label, outputPosition=NODAL)
viewportObj.odbDisplay.contourOptions.setValues(maxAutoCompute=ON, minAutoCompute=ON)
full_cartilage(label)

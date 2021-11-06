import numpy as np
import pickle
import tensorflow as tf
import pyvista as pv

# %% Enable GPU Memory Growth
def enableGPUMemGro():
    gpus = tf.config.experimental.list_physical_devices('GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
        
# %% Helper Functions for Data Loader

# Ususal tf dtype
fdtype = tf.float32

# Ones like mesh shape
point_to_field = lambda x, shape: x * np.ones([*shape,1])

# A pickle file loader from a tf tensor
pklLoader = lambda x: pickle.load(open(x.numpy(), 'rb'))

# Concatenator for a tf tensors
concatenator = lambda a,b: np.concatenate((a.numpy(),b.numpy()), axis=-1)

# Standardize
standardizer = lambda x, mean, std: (x.numpy() - mean) / (std + 1e-16)

# Relative Error in L1 norm
def relErrL1(y_true , y_pred):
    SMALL = 1e-16
    y_true, y_pred = tf.cast(y_true, fdtype), tf.cast(y_pred, fdtype)
    relErr = tf.divide(tf.abs(y_true-y_pred), tf.abs(y_true) + SMALL)
    return tf.reduce_mean(relErr)

# Relative Error in L2 norm
def relErrL2(y_true , y_pred):
    SMALL = 1e-16
    y_true, y_pred = tf.cast(y_true, fdtype), tf.cast(y_pred, fdtype)
    relErr = tf.divide(tf.norm(y_true-y_pred), tf.norm(y_true) + SMALL)
    return relErr

# %% Data Loader

# U_mean, U_std = 0.,1.
# TI_mean, TI_std = 0.,1. 

U_mean, U_std = 7.3, 2.3
TI_mean, TI_std = 10.3, 4.8

getU = lambda defU, UHub: UHub*(1-defU)
getTI = lambda TIAdd, TIHub: TIAdd + TIHub
    
# Hub Data (Field)
def loadPickledSamplesAtInputFields(s, meshShape):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, U_mean, U_std], [fdtype])
    UHub_field = tf.py_function(point_to_field, [UHub, meshShape], [fdtype])  
    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TI_mean, TI_std], [fdtype])
    TIHub_field = tf.py_function(point_to_field, [TIHub, meshShape], [fdtype])
    
    data = tf.py_function(concatenator, [UHub_field, TIHub_field], [fdtype])
    data = tf.reshape(data, [*meshShape,2])
    # return data
    
    # AData = tf.py_function(pklLoader, [s+'/A.pkl'], [fdtype])
    AData = tf.py_function(pklLoader, [s+'/R.pkl'], [fdtype])
    AData = tf.reshape(AData, [*meshShape,6])
    # return AData
    
    conc_data = tf.py_function(concatenator, [AData, data], [fdtype])
    conc_data = tf.reshape(conc_data, [*meshShape,8])
    return conc_data

# Hub Data (Point)
def loadPickledSamplesAtHub(s):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, U_mean, U_std], [fdtype])    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TI_mean, TI_std], [fdtype])
    data = tf.py_function(concatenator, [UHub, TIHub], [fdtype])
    data = tf.reshape(data, [2])
    return data

# Anisotropy Data
def loadPklSamplesAData(s, meshShape):
    data = tf.py_function(pklLoader, [s+'/A.pkl'], [fdtype])
    data = tf.reshape(data, [*meshShape,6])
    return data

# Output Data
def loadPklSamplesOutputFields(s, meshShape):
    defU = tf.py_function(pklLoader, [s+'/defU.pkl'], [fdtype])
    TIAdd = tf.py_function(pklLoader, [s+'/TIAdd.pkl'], [fdtype])
    
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    
    UMag = tf.py_function(getU, [defU, UHub], [fdtype])
    TI = tf.py_function(getTI, [TIAdd, TIHub], [fdtype])
    
    data = tf.py_function(concatenator, [UMag, TI], [fdtype])
    # data = tf.py_function(concatenator, [defU, TIAdd], [fdtype])
    data = tf.reshape(data, [*meshShape,2])
    return data

# Load Data
def loadData(fileList, meshShape):
    hubData = fileList.map(lambda x: loadPickledSamplesAtHub(x)) 
    input_fields = fileList.map(
        lambda x: loadPickledSamplesAtInputFields(x, meshShape)
    )
    AData = fileList.map(lambda x: loadPklSamplesAData(x, meshShape))
    output_fields = fileList.map(lambda x: loadPklSamplesOutputFields(x, meshShape))
    return hubData, AData, input_fields, output_fields

# Split and Batch Data
def batchSplitData(data, trainFrac, batchSize):
    if trainFrac < 1.:
        train_size = int(len(data)*trainFrac)
        train_data = data.take(train_size)
        rest_data  = data.skip(train_size)
    
        valid_size = int(len(rest_data)/2.)
        valid_data = rest_data.take(valid_size)
        test_data  = rest_data.skip(valid_size)
    
        train_data_batched = train_data.batch(batchSize).cache().prefetch(1)
        valid_data_batched = valid_data.batch(batchSize).cache().prefetch(1)
        # test_data_batched = test_data.batch(batchSize).cache().prefetch(1)
        test_data_batched = test_data.batch(1).cache().prefetch(1)
    
        return train_data_batched, valid_data_batched, test_data_batched
    
    return data.batch(1).cache().prefetch(1)

# Data Generator Class
class dataGenerator():
    def __init__(self, dataDir, meshName, sampleRegx,
                 meshShape, trainFrac=0.8, batchSize=1):
        
        self.batchSize = batchSize
        self.meshShape = meshShape
        self.trainFrac = trainFrac
        self.dataDir = dataDir+meshName+sampleRegx
        
        self.fileList = tf.data.Dataset.list_files(self.dataDir, shuffle=False)
        self.hubData, self.AData, self.input_fields, self.output_fields = \
            loadData(self.fileList, self.meshShape)
        
        # Dat for UNet Model
        self.UNetIOData = tf.data.Dataset.zip(
            (self.input_fields, self.output_fields)
        )
        self.UNetIOBatchedSplitData = batchSplitData(
            self.UNetIOData, self.trainFrac, self.batchSize
        )
        
        # Data for UNetAug Model
        input_ = tf.data.Dataset.zip((self.AData, self.hubData))
        self.UNetAugIOData = tf.data.Dataset.zip(
            (input_, self.output_fields)
        )
        self.UNetAugIOBatchedSplitData = batchSplitData(
            self.UNetAugIOData, self.trainFrac, self.batchSize
        )

# %%
import os
import numpy as np
import pickle
import tensorflow as tf
import random

# %% Global Seeding for Reproducibility
def set_global_determinism(seed=42, fast_n_close=False):
    """
        Enable 100% reproducibility on operations related to 
        tensor and randomness.
        https://suneeta-mall.github.io/2019/12/22/Reproducible-ml-tensorflow.html

        Parameters:
        seed (int): seed value for global randomness
        fast_n_close (bool): whether to achieve efficient at 
        the cost of determinism/reproducibility
    """
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    tf.random.set_seed(seed)
    np.random.seed(seed)
    
    if fast_n_close:
        return
    
    os.environ['TF_DETERMINISTIC_OPS'] = '1'
    os.environ['TF_CUDNN_DETERMINISTIC'] = '1'
    tf.config.threading.set_inter_op_parallelism_threads(1)
    tf.config.threading.set_intra_op_parallelism_threads(1)

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

# Relative Error in max norm
def max_norm(y_true , y_pred):
    SMALL = 1e-16
    y_true, y_pred = tf.cast(y_true, fdtype), tf.cast(y_pred, fdtype)
    relErr = tf.divide(tf.abs(y_true-y_pred), tf.abs(y_true) + SMALL)
    return relErr.max()

# Relative Error in L1 norm
def L1(y_true , y_pred):
    SMALL = 1e-16
    y_true, y_pred = tf.cast(y_true, fdtype), tf.cast(y_pred, fdtype)
    relErr = tf.divide(tf.abs(y_true-y_pred), tf.abs(y_true) + SMALL)
    return tf.reduce_mean(relErr)

# Relative Error in L2 norm
def L2(y_true , y_pred):
    SMALL = 1e-16
    y_true, y_pred = tf.cast(y_true, fdtype), tf.cast(y_pred, fdtype)
    relErr = tf.divide(tf.norm(y_true-y_pred), tf.norm(y_true) + SMALL)
    return relErr

# %% Data Loader

# From the Data Processing Step
UHub_mean, UHub_std = (6.279, 1.967)
TIHub_mean, TIHub_std = (12.969, 4.438)
UMagMax, TIMax, tkeMax = (14.616, 21.795, 8.815)

# Hub Data (Field)
def loadPickledSamplesAtInputFields(s, meshShape):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, UHub_mean, UHub_std], [fdtype])
    UHub_field = tf.py_function(point_to_field, [UHub, meshShape], [fdtype])  
    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHub_mean, TIHub_std], [fdtype])
    TIHub_field = tf.py_function(point_to_field, [TIHub, meshShape], [fdtype])
    
    dataHub = tf.py_function(concatenator, [UHub_field, TIHub_field], [fdtype])
    dataHub = tf.reshape(dataHub, [*meshShape,2])
    
    # AData = tf.py_function(pklLoader, [s+'/A.pkl'], [fdtype])
    # AData = tf.reshape(AData, [*meshShape,6])
    
    # conc_data = tf.py_function(concatenator, [AData, dataHub], [fdtype])
    # conc_data = tf.reshape(conc_data, [*meshShape,8])
    
    C_vec = tf.py_function(pklLoader, [s+'/C_vec.pkl'], [fdtype])
    C_vec = tf.reshape(C_vec, [*meshShape,2])
    
    conc_data = tf.py_function(concatenator, [C_vec, dataHub], [fdtype])
    conc_data = tf.reshape(conc_data, [*meshShape,4])

    return conc_data

# Hub Data (Point)
def loadPickledSamplesAtHub(s):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, UHub_mean, UHub_std], [fdtype])    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHub_mean, TIHub_std], [fdtype])
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
    UMag = tf.py_function(pklLoader, [s+'/UMag.pkl'], [fdtype])
    TI = tf.py_function(pklLoader, [s+'/TI.pkl'], [fdtype])
    data = tf.py_function(concatenator, [UMag, TI], [fdtype])
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
        test_data_batched = test_data.batch(1).cache().prefetch(1)
    
        return train_data_batched, valid_data_batched, test_data_batched
    
    return data.batch(1).cache().prefetch(1)

# Data Generator Class
class dataGenerator():
    def __init__(self, fileNames, meshShape, trainFrac=0.8, batchSize=2):
        
        # Attributes
        self.batchSize = batchSize
        self.meshShape = meshShape
        self.trainFrac = trainFrac
        
        # Load Files and Data
        self.fileList = tf.data.Dataset.from_tensor_slices(fileNames)
        self.hubData, self.AData, self.input_fields, self.output_fields = \
            loadData(self.fileList, self.meshShape)
        
        # Data for UNet Model
        self.UNetIOData = tf.data.Dataset.zip(
            (self.input_fields, self.output_fields)
        )
        self.UNetIOBatchedSplitData = batchSplitData(
            self.UNetIOData, self.trainFrac, self.batchSize
        )
        
import os
import numpy as np
import pickle
import tensorflow as tf
import random

# %% Global Seeding for Reproducibility
def set_global_determinism(seed=42, fast_n_close=False):
    """
        Enable 100% reproducibility on operations related to 
        tensor and randomness. Reproducibility is gauteed only
        when the entire script is ran again from start.
        
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
fdtype = tf.float64

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

# Defined in dataGenerator class
UHub_mean, UHub_std = (None, None)
TIHub_mean, TIHub_std = (None, None)

# Hub Data (Point data at hub height)
def loadUHubTIHubSamples(s):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, UHub_mean, UHub_std], [fdtype])    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHub_mean, TIHub_std], [fdtype])
    data = tf.py_function(concatenator, [UHub, TIHub], [fdtype])
    data = tf.reshape(data, [2])
    return data

# Anisotropy Data (Aij or C_vec)
def loadAnisoSamples(s, meshShape):
    data = tf.py_function(pklLoader, [s+'/A.pkl'], [fdtype])
    data = tf.reshape(data, [*meshShape,6])
    
    # data = tf.py_function(pklLoader, [s+'/C_vec.pkl'], [fdtype])
    # data = tf.reshape(data, [*meshShape,2])
    
    return data

# Concat [ Hub (Point Data -> Field) | Anisotropy Data (Aij or C_vec) ]
def loadAllInputFieldSamples(s, meshShape):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, UHub_mean, UHub_std], [fdtype])
    UHub_field = tf.py_function(point_to_field, [UHub, meshShape], [fdtype])  
    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHub_mean, TIHub_std], [fdtype])
    TIHub_field = tf.py_function(point_to_field, [TIHub, meshShape], [fdtype])
    
    dataHub = tf.py_function(concatenator, [UHub_field, TIHub_field], [fdtype])
    dataHub = tf.reshape(dataHub, [*meshShape,2])
    
    AData = tf.py_function(pklLoader, [s+'/A.pkl'], [fdtype])
    AData = tf.reshape(AData, [*meshShape,6])
    
    conc_data = tf.py_function(concatenator, [AData, dataHub], [fdtype])
    conc_data = tf.reshape(conc_data, [*meshShape,8])
    
    # C_vec = tf.py_function(pklLoader, [s+'/C_vec.pkl'], [fdtype])
    # C_vec = tf.reshape(C_vec, [*meshShape,2])
    
    # conc_data = tf.py_function(concatenator, [C_vec, dataHub], [fdtype])
    # conc_data = tf.reshape(conc_data, [*meshShape,4])

    return conc_data

# Output Data
def loadAllOutputFieldSamples(s, meshShape):
    UMag = tf.py_function(pklLoader, [s+'/UMag.pkl'], [fdtype])
    TI = tf.py_function(pklLoader, [s+'/TI.pkl'], [fdtype])
    data = tf.py_function(concatenator, [UMag, TI], [fdtype])
    data = tf.reshape(data, [*meshShape,2])
    return data

# Load Data
def loadData(fileList, meshShape):
    input_fields = fileList.map(
        lambda x: loadAllInputFieldSamples(x, meshShape)
    )
    output_fields = fileList.map(
        lambda x: loadAllOutputFieldSamples(x, meshShape)
    )
    UHubTIHubData = fileList.map(lambda x: loadUHubTIHubSamples(x)) 
    anisoData = fileList.map(lambda x: loadAnisoSamples(x, meshShape))
    return UHubTIHubData, anisoData, input_fields, output_fields

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
    def __init__(self, mlDataDir, mlMeshName, fileNames, 
                 meshShape, trainFrac=0.8, batchSize=2):
        
        # From the Data Processing Step
        global UHub_mean, UHub_std, TIHub_mean, TIHub_std
        statsDir = mlDataDir+mlMeshName
        UHub_mean = pickle.load(open(statsDir+'UHubMean.pkl','rb'))
        UHub_std = pickle.load(open(statsDir+'UHubStd.pkl','rb'))
        TIHub_mean = pickle.load(open(statsDir+'TIHubMean.pkl','rb'))
        TIHub_std = pickle.load(open(statsDir+'TIHubStd.pkl','rb'))
        
        # Attributes
        self.batchSize, self.meshShape, self.trainFrac = \
            batchSize, meshShape, trainFrac
        self.UHub_mean, self.UHub_std, self.TIHub_mean, self.TIHub_std = \
            UHub_mean, UHub_std, TIHub_mean, TIHub_std
        
        # Load Files and Data 
        self.fileList = tf.data.Dataset.from_tensor_slices(fileNames)
        self.UHubTIHubData, self.anisoData, \
            self.input_fields, self.output_fields = \
                loadData(self.fileList, self.meshShape)
        
        # Data for UNet Model: zip(concat(anisoField,hubFields),outFields)
        self.UNetIOData = tf.data.Dataset.zip(
            (self.input_fields, self.output_fields)
        )
        self.UNetIOBatchedSplitData = batchSplitData(
            self.UNetIOData, self.trainFrac, self.batchSize
        )
        
        # Data for AutoEnc Model: anisoData
        self.AutoEncIOData = tf.data.Dataset.zip(
            (self.anisoData, self.anisoData)
        )
        self.AutoEncIOBatchedSplitData = batchSplitData(
            self.AutoEncIOData, self.trainFrac, self.batchSize
        )
        
        # Data for TransposedCNN Model: zip(zip(anisoField,hubData),outFields)
        self.TransposedCNNIOData = tf.data.Dataset.zip(
            ( 
                tf.data.Dataset.zip((self.anisoData, self.UHubTIHubData)),
                self.output_fields
            )
        )
        self.TransposedCNNIOBatchedSplitData = batchSplitData(
            self.TransposedCNNIOData, self.trainFrac, self.batchSize
        )
        
import numpy as np
import pickle
import tensorflow as tf

# %% Enable GPU Memory Growth
def enableGPUMemGro():
    gpus = tf.config.experimental.list_physical_devices('GPU')
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)

# %% Helper Functions

# Ususal tf dtype
fdtype = tf.float32

# A pickle file loader from a tf tensor
pklLoader = lambda x: pickle.load(open(x.numpy(), 'rb'))

# Concatenator for a tf tensors
concatenator = lambda a,b: np.concatenate((a.numpy(),b.numpy()), axis=-1)

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

# Hub Data
def loadPickledSamplesAtHub(s):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
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
    data = tf.py_function(concatenator, [defU, TIAdd], [fdtype])
    data = tf.reshape(data, [*meshShape,2])
    return data

# Load Data
def loadData(fileList, meshShape):
    hubData = fileList.map(
        lambda x: loadPickledSamplesAtHub(x), deterministic=False)

    AData = fileList.map(
        lambda x: loadPklSamplesAData(x,meshShape), deterministic=False)

    output = fileList.map(
        lambda x: loadPklSamplesOutputFields(x,meshShape), deterministic=False)

    input_ = tf.data.Dataset.zip((AData, hubData))
    return tf.data.Dataset.zip((input_, output))

# Split and Batch Data
def batchSplitData(ioData, trainFrac, batchSize):
    train_size = int(len(ioData)*trainFrac)
    train_data = ioData.take(train_size)
    rest_data  = ioData.skip(train_size)

    valid_size = int(len(rest_data)/2.)
    valid_data = rest_data.take(valid_size)
    test_data  = rest_data.skip(valid_size)

    train_data_batched = train_data.batch(batchSize).cache().prefetch(1)
    valid_data_batched = valid_data.batch(batchSize).cache().prefetch(1)
    test_data_batched = test_data.batch(batchSize).cache().prefetch(1)

    return train_data_batched, valid_data_batched, test_data_batched

# Data Generator Class
class dataGenerator():
    def __init__(self, dataDir, meshName, sampleRegx,
                 meshShape, trainFrac=0.8, batchSize=1):
        self.batchSize = batchSize
        self.meshShape = meshShape
        self.trainFrac = trainFrac
        self.dataDir = dataDir+meshName+sampleRegx
        self.fileList = tf.data.Dataset.list_files(self.dataDir, shuffle=False)
        self.ioData = loadData(self.fileList, self.meshShape)
        self.ioBatchedSplitData = batchSplitData(
            self.ioData, self.trainFrac, self.batchSize)

# %%
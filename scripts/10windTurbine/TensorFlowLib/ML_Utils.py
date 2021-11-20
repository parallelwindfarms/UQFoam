import os
import numpy as np
import pickle
import tensorflow as tf
import random
import matplotlib.pyplot as plt

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
        
# %% Plotting True vs Pred vs |True-Pred|
def makePlots(s, mlMeshShape, y0Plane_WT_idx, zhPlane_WT_idx, 
              UMagTestTrue, UMagTestPred, UMagDiff,
              TITestTrue, TITestPred, TIDiff):
    plot_soln = lambda ax, soln, plane, norm=None: ax.imshow(
        soln[plane].reshape(-1,mlMeshShape[2])[::-1], 
        aspect='auto', 
        norm=norm
    )
    err_pct = 0.05
    normGlobal = plt.Normalize(0, 0.05)

    # Contour Plots
    fig, ax = plt.subplots(ncols=3, nrows=4, constrained_layout=True, 
                           sharex=True, figsize=(16,6))
    ax, CS = ax.flat, [0]*12
    
    norm = plt.Normalize(UMagTestTrue[y0Plane_WT_idx].min(), UMagTestTrue[y0Plane_WT_idx].max())
    CS[0] = plot_soln(ax[0], UMagTestTrue, y0Plane_WT_idx)
    CS[1] = plot_soln(ax[1], UMagTestPred, y0Plane_WT_idx, norm=norm)
    norm = plt.Normalize(0, UMagTestTrue[y0Plane_WT_idx].max()*err_pct)
    CS[2] = plot_soln(ax[2], UMagDiff, y0Plane_WT_idx, norm=normGlobal)
    
    CS[3] = plot_soln(ax[3], TITestTrue, y0Plane_WT_idx)
    CS[4] = plot_soln(ax[4], TITestPred, y0Plane_WT_idx)
    norm = plt.Normalize(0, TITestTrue[y0Plane_WT_idx].max()*err_pct)
    CS[5] = plot_soln(ax[5], TIDiff, y0Plane_WT_idx, norm=normGlobal)
    
    norm = plt.Normalize(UMagTestTrue[zhPlane_WT_idx].min(), UMagTestTrue[zhPlane_WT_idx].max())
    CS[6] = plot_soln(ax[6], UMagTestTrue, zhPlane_WT_idx)
    CS[7] = plot_soln(ax[7], UMagTestPred, zhPlane_WT_idx, norm=norm)
    norm = plt.Normalize(0, UMagTestTrue[zhPlane_WT_idx].max()*err_pct)
    CS[8] = plot_soln(ax[8], UMagDiff, zhPlane_WT_idx, norm=normGlobal)
    
    CS[9] = plot_soln(ax[9], TITestTrue, zhPlane_WT_idx)
    CS[10] = plot_soln(ax[10], TITestPred, zhPlane_WT_idx)
    norm = plt.Normalize(0, TITestTrue[zhPlane_WT_idx].max()*err_pct)
    CS[11] = plot_soln(ax[11], TIDiff, zhPlane_WT_idx, norm=normGlobal)
        
    for i in range(12):
        # ax[i].set_xticks([])
        ax[i].set_yticks([])
        if i in [1,4,7,10]:
            fig.colorbar(CS[i-1], ax=ax[i], aspect=50)
        else:
            fig.colorbar(CS[i], ax=ax[i], aspect=50)
    
    ax[0].set_title('OpenFOAM (True)')
    ax[1].set_title('U-Net (Pred)')
    ax[2].set_title('|True - Pred|')
    ax[0].set_ylabel('UMag')
    ax[3].set_ylabel('TI')
    ax[6].set_ylabel('UMag')
    ax[9].set_ylabel('TI')
    
    # ax[0].set_xlim(64*0,64*3)
    
    fig.suptitle(f'Case #{s}')
        
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
UHubMean, UHubStd = (None, None)
TIHubMean, TIHubStd = (None, None)

# Hub Data (Point data at hub height)
def loadUHubTIHubSamples(s):
    UHub = tf.py_function(pklLoader, [s+'/UHub.pkl'], [fdtype])
    UHub = tf.py_function(standardizer, [UHub, UHubMean, UHubStd], [fdtype])    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHubMean, TIHubStd], [fdtype])
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
    UHub = tf.py_function(standardizer, [UHub, UHubMean, UHubStd], [fdtype])
    UHub_field = tf.py_function(point_to_field, [UHub, meshShape], [fdtype])  
    
    TIHub = tf.py_function(pklLoader, [s+'/TIHub.pkl'], [fdtype])
    TIHub = tf.py_function(standardizer, [TIHub, TIHubMean, TIHubStd], [fdtype])
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
        global UHubMean, UHubStd, TIHubMean, TIHubStd
        statsDir = mlDataDir+mlMeshName
        UHubMean = pickle.load(open(statsDir+'UHubMean.pkl','rb'))
        UHubStd = pickle.load(open(statsDir+'UHubStd.pkl','rb'))
        TIHubMean = pickle.load(open(statsDir+'TIHubMean.pkl','rb'))
        TIHubStd = pickle.load(open(statsDir+'TIHubStd.pkl','rb'))
        
        # Attributes
        self.batchSize, self.meshShape, self.trainFrac = \
            batchSize, meshShape, trainFrac
        self.UHubMean, self.UHubStd = UHubMean, UHubStd
        self.TIHubMean, self.TIHubStd = TIHubMean, TIHubStd
        
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
        
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Input, Conv3D, Conv3DTranspose
from tensorflow.keras.layers import concatenate, UpSampling3D, Reshape
from tensorflow.keras.layers import Dropout, BatchNormalization
from tensorflow.keras.layers import ReLU, LeakyReLU, ELU
from tensorflow.keras.regularizers import l1

# %% Conv Block
def convBlock(name='NO_NAME', nOutChannels=None, convType=None, kernelSize=4,
              stride=2, activation='relu', kernelInitializer='he_normal',
              padding='same', batchNorm=True, dropFrac=0., l1_lambda=0.):

    if activation=='relu':
        activation = ReLU()
    elif activation=='leaky_relu':
        activation = LeakyReLU()
    elif activation=='elu':
        activation = ELU()
    else:
        activation = None

    block = Sequential(name=name)

    if convType=='upsampled':
        block.add(UpSampling3D(size=2))
        block.add(Conv3D(nOutChannels, kernelSize-1, strides=1,
                         activation=activation, padding=padding,
                         kernel_initializer=kernelInitializer,
                         kernel_regularizer=l1(l1_lambda)))
    elif convType=='transposed':
        block.add(Conv3DTranspose(nOutChannels, kernelSize, strides=stride,
                                  activation=activation, padding=padding,
                                  kernel_initializer=kernelInitializer))
    else:
        block.add(Conv3D(nOutChannels, kernelSize, strides=stride,
                         activation=activation, padding=padding,
                         kernel_initializer=kernelInitializer))

    if batchNorm:
        block.add(BatchNormalization())
    if dropFrac>0.:
        block.add(Dropout(dropFrac))

    return block

# %% AutoEnc Model
def AutoEnc(meshShape, nRandFieldsChannels=2, nOutChannels=2, 
            nLatentChannels=10, dropFrac=0., channels=64, 
            l1_lambda=0., convType='upsampled', latent=True):

    # Clear Session
    tf.keras.backend.clear_session()

    # Check for a valid mesh shape
    assert (meshShape[0]==64 or meshShape[0]==128)

    # Params
    act = 'leaky_relu'

    # Inputs Layers
    randFields = Input(shape=(*meshShape, nRandFieldsChannels), name='input')

    # Hidden Layers
    if meshShape[0]==128:
        cLayer_i = convBlock('cLayer_i', channels, activation=None, dropFrac=0., batchNorm=False, l1_lambda=l1_lambda)
        cLayer_0 = convBlock('cLayer_0', channels*2, activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)
    elif meshShape[0]==64:
        cLayer_0 = convBlock('cLayer_0', channels*2, activation=None, dropFrac=0., batchNorm=False, l1_lambda=l1_lambda)

    cLayer_1 = convBlock('cLayer_1', channels*2, activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)
    cLayer_2 = convBlock('cLayer_2', channels*4, activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)
    cLayer_3 = convBlock('cLayer_3', channels*8, activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)
    cLayer_4 = convBlock('cLayer_4', channels*8, activation=act, dropFrac=dropFrac, kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    cLayer_5 = convBlock('cLayer_5', channels*8, activation=act, dropFrac=dropFrac, kernelSize=2, padding='valid', batchNorm=False, l1_lambda=l1_lambda)

    dLayer_5 = convBlock('dLayer_5', channels*8, convType=convType, dropFrac=dropFrac, kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    dLayer_4 = convBlock('dLayer_4', channels*8, convType=convType, dropFrac=dropFrac, kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    dLayer_3 = convBlock('dLayer_3', channels*4, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)
    dLayer_2 = convBlock('dLayer_2', channels*2, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)
    dLayer_1 = convBlock('dLayer_1', channels*2, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)

    if meshShape[0]==128:
        dLayer_0 = convBlock('dLayer_0', channels, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)

    # Output Layer
    dLayer_o = convBlock('dLayer_o', nOutChannels, activation=None, convType='transposed', dropFrac=0., batchNorm=False, l1_lambda=l1_lambda)

    # Model
    if meshShape[0]==128:
        cOut_i = cLayer_i(randFields)
        cOut_0 = cLayer_0(cOut_i)
    elif meshShape[0]==64:
        cOut_0 = cLayer_0(randFields)

    cOut_1 = cLayer_1(cOut_0)
    cOut_2 = cLayer_2(cOut_1)
    cOut_3 = cLayer_3(cOut_2)
    cOut_4 = cLayer_4(cOut_3)
    cOut_5 = cLayer_5(cOut_4)
    
    latentSpace = convBlock('latentSpace', nLatentChannels, activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)(cOut_5)
    latentOut = convBlock('latentOut', cOut_5.shape[-1], activation=act, dropFrac=dropFrac, l1_lambda=l1_lambda)(latentSpace)

    if latent:
        dOut_5_cOut_4 = concatenate([dLayer_5(latentOut), cOut_4], name='dLayer_5_conc')
    else:
        dOut_5_cOut_4 = concatenate([dLayer_5(cOut_5), cOut_4], name='dLayer_5_conc')
    
    dOut_4_cOut_3 = concatenate([dLayer_4(dOut_5_cOut_4), cOut_3], name='dLayer_4_conc')
    dOut_3_cOut_2 = concatenate([dLayer_3(dOut_4_cOut_3), cOut_2], name='dLayer_3_conc')
    dOut_2_cOut_1 = concatenate([dLayer_2(dOut_3_cOut_2), cOut_1], name='dLayer_2_conc')
    dOut_1_cOut_i = concatenate([dLayer_1(dOut_2_cOut_1), cOut_0], name='dLayer_1_conc')

    if meshShape[0]==128:
        dOut_0_cOut_i = concatenate([dLayer_0(dOut_1_cOut_i), cOut_i], name='dLayer_0_conc')
        dOut_o  = dLayer_o(dOut_0_cOut_i)
    elif meshShape[0]==64:
        dOut_o = dLayer_o(dOut_1_cOut_i)

    model = tf.keras.Model(inputs=[randFields], outputs=[dOut_o], name='AutoEnc')

    return model

# %% Latent Space Representation
def LatentSpaceRepr(model):
    latentVars = Sequential()
    for layer in model.layers:
        latentVars.add(layer)
        if layer.name == 'latentSpace':
            break
    latentVars.trainable = False
    return latentVars

# %% TransposedCNN Model
def TransposedCNN(meshShape, autoencoder, nRandFields=2, nHubData=2, 
                  nOutChannels=2, dropFrac=0., channels=64, l1_lambda=0., 
                  convType='transposed'):

    # Clear Session
    tf.keras.backend.clear_session()

    # Check for a valid mesh shape
    assert (meshShape[0]==64 or meshShape[0]==128)

    # Defining Inputs
    randomFields = Input(shape=(*meshShape, nRandFields), name='randomFields')
    hubData = Input(shape=(nHubData), name='hubData')
    
    # Inputs Layer
    latentSpaceShape = autoencoder.get_layer('latentSpace').output.shape[1:]
    latentVars = LatentSpaceRepr(autoencoder)(randomFields)
    # hubDataReshaped = Reshape((*latentSpaceShape[:-1],-1))(hubData)
    
    hubData_ = Reshape((*latentSpaceShape[:-1],-1))(hubData)
    hubDataReshaped = hubData_
    for i in range(latentSpaceShape[-1]-1):
        hubDataReshaped = concatenate([hubDataReshaped, hubData_])
        
    inputLayer = concatenate([latentVars, hubDataReshaped], name='inputLayer')
    
    # Hidden Layers
    dLayer_7 = convBlock('dLayer_7', channels*8, activation=None, dropFrac=0., batchNorm=False, l1_lambda=0.)
    dLayer_6 = convBlock('dLayer_6', channels*16, convType=convType, dropFrac=dropFrac,  kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    dLayer_5 = convBlock('dLayer_5', channels*16, convType=convType, dropFrac=dropFrac, kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    dLayer_4 = convBlock('dLayer_4', channels*8, convType=convType, dropFrac=dropFrac, kernelSize=2, padding='valid', l1_lambda=l1_lambda)
    dLayer_3 = convBlock('dLayer_3', channels*4, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)
    dLayer_2 = convBlock('dLayer_2', channels*4, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)
    if meshShape[0]==128:
        dLayer_1 = convBlock('dLayer_0', channels, convType=convType, dropFrac=dropFrac, l1_lambda=l1_lambda)
        
    # Output Layer
    dLayer_o = convBlock('dLayer_o', nOutChannels, convType='transposed', dropFrac=0., batchNorm=False, l1_lambda=l1_lambda)

    # Model
    dOut_7 = dLayer_7(inputLayer)
    dOut_6 = dLayer_6(dOut_7)
    dOut_5 = dLayer_5(dOut_6)
    dOut_4 = dLayer_4(dOut_5)
    dOut_3 = dLayer_3(dOut_4)
    dOut_2 = dLayer_2(dOut_3)
    
    if meshShape[0]==128:
        dOut_0 = dLayer_1(dOut_2)
        dOut_o = dLayer_o(dOut_0)
    elif meshShape[0]==64:
        dOut_o = dLayer_o(dOut_2)
    
    model = tf.keras.Model(inputs=[randomFields, hubData], outputs=[dOut_o], name='TransposedCNN')

    return model
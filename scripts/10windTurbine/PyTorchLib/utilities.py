import numpy as np
import torch
import torch.nn as nn

import sys
from subprocess import call

#%% Cuda info
def cudaInfo():
    print('__Python VERSION:', sys.version)
    print('__pyTorch VERSION:', torch.__version__)
    print('__CUDA VERSION')
    print('__CUDNN VERSION:', torch.backends.cudnn.version())
    print('__Number CUDA Devices:', torch.cuda.device_count())
    print('__Devices')
    call(["nvidia-smi", "--format=csv", "--query-gpu=index,name,driver_version,memory.total,memory.used,memory.free"])
    print('Active CUDA Device: GPU', torch.cuda.current_device())

    print ('Available devices ', torch.cuda.device_count())
    print ('Current cuda device ', torch.cuda.current_device())

#%% grad, jacobian, divergence, laplace
def grad(f, z, retain_graph=True, create_graph=True, allow_unused=False):
    return torch.autograd.grad(f, z, retain_graph=retain_graph,
                               create_graph=create_graph,
                               allow_unused=allow_unused,
                               grad_outputs=torch.ones_like(f))

def gradient(y, x, grad_outputs=None):
    """Compute dy/dx @ grad_outputs"""
    if grad_outputs is None:
        grad_outputs = torch.ones_like(y)
    grad = torch.autograd.grad(y, [x], grad_outputs = grad_outputs, create_graph=True)[0]
    return grad

def jacobian(y, x):
    """Compute dy/dx = dy/dx @ grad_outputs;
    for grad_outputs in [1, 0, ..., 0], [0, 1, 0, ..., 0], ...., [0, ..., 0, 1]"""
    jac = torch.zeros(y.shape[0], x.shape[0])
    for i in range(y.shape[0]):
        grad_outputs = torch.zeros_like(y)
        grad_outputs[i] = 1
        jac[i] = gradient(y, x, grad_outputs = grad_outputs)
    return jac

def divergence(y, x):
    div = 0.
    for i in range(y.shape[-1]):
        div += torch.autograd.grad(y[..., i], x, torch.ones_like(y[..., i]), create_graph=True)[0][..., i:i+1]
    return div

def laplace(y, x):
    grad = gradient(y, x)
    return divergence(grad, x)

#%% L2loss
def L2Loss(net, device):
    loss = torch.tensor(0.).to(device)
    for m in net.modules():
        if hasattr(m, 'weight'):
            loss += m.weight.norm()**2
    return loss

#%% Total parameters
def totalParams(net):
    s = sum(p.numel() for p in net.parameters() if p.requires_grad)
    print(f'Number of trainable parameters = {s} ')

#%% splitIndices
def splitIndices(n, valPct):
    # size of val set
    nVal = int(valPct*n)
    # random permutations of 0 to n-1
    idxs = np.random.permutation(n)
    # pick first n_val indices for val set
    return idxs[nVal:], idxs[:nVal]

#%% Get CPU/GPU device
def getDefaultDevice(gpu=0, cpu=False):
    """Pick GPU if available, else CPU"""
    if torch.cuda.is_available():
        return torch.device('cuda:'+str(gpu))
    return torch.device('cpu')

#%% Move model or data to chosen device
def toDevice(data, device):
    """Move tensor(s) to chosen device"""
    if isinstance(data, (list,tuple)):
        return [toDevice(x, device) for x in data]
    return data.to(device, non_blocking=True)

#%% Move data to device on the fly (only when it is about to be used)
class DeviceDataLoader():
    """Wrap a dataloader to move data to a device"""
    def __init__(self, dl, device):
        self.dl = dl
        self.device = device

    def __iter__(self):  # called when dl is used in a for loop
        """Yield a batch of data after moving it to device"""
        for b in self.dl:
            yield toDevice(b, self.device) # return when needed (called)

    def __len__(self):
        """Number of batches"""
        return len(self.dl)

#%% weightsInitXavier
def weightsInitXavier(m):
    if type(m) == nn.Linear:
        torch.nn.init.xavier_uniform_(m.weight)
        m.bias.data.fill_(0.01)

#%% weightsInitUniform
def weightsInitUniform(m):
    classname = m.__class__.__name__
    # for every Linear layer in a model..
    if classname.find('Linear') != -1:
        # get the number of the inputs
        n = m.in_features
        y = 1.0/np.sqrt(n)
        m.weight.data.uniform_(-y, y)
        m.bias.data.fill_(0)

#%% randomToDevice
def randomToDevice(nrows=1, ncols=1, device='gpu'):
    return (torch.rand((nrows, ncols), requires_grad=True)).to(device)

#%% Network w/ output
def Net(net, X, nOutCells):
    soln = net(X)
    if len(soln.shape)>=2:
        UMag = [soln[s][nOutCells*0:nOutCells*1] for s in range(soln.shape[0])]
        tke  = [soln[s][nOutCells*1:nOutCells*2] for s in range(soln.shape[0])]
        return UMag, tke
    else:
        UMag = soln[nOutCells*0:nOutCells*1]
        tke  = soln[nOutCells*1:nOutCells*2]
        return UMag, tke
    
#%% Torch L1 Error
def L1Err(y_true, y_pred):
    SMALL = 1e-16    
    err = torch.abs(y_true-y_pred)/(torch.abs(y_true)+SMALL)
    return round(err.mean().item(), 3)
    # return round(err.max().item(), 3)

#%% Torch L2 Error
def L2Err(y_true, y_pred):
    SMALL = 1e-16    
    err = torch.norm(y_true-y_pred)/(torch.norm(y_true)+SMALL)
    return round(err.item(), 3)

#%% To numpy
torchToNumpy = lambda x: x.detach().cpu()#.numpy()
    






import torch
import torch.nn as nn

#%% A template FCNN in Pytorch
class NeuralNetwork(nn.Module):
    """
    An implementation of a fully connected feed forward
    Neural network in pytorch.
    """
    def __init__(self, layersizes=[1, 1],
                 activation=torch.relu,
                 final_layer_activation=None):
        super(NeuralNetwork, self).__init__()
        self.layersizes = layersizes
        self.input_dim = self.layersizes[0]
        self.hidden_sizes = self.layersizes[1:-1]
        self.output_dim = self.layersizes[-1]
        self.activation = activation
        self.final_layer_activation = final_layer_activation
        if self.final_layer_activation is None:
            self.final_layer_activation = nn.Identity()
        self.nlayers = len(self.hidden_sizes) + 1
        self.layernames = [] ## Dictionary to store all the FC layers

        # define FC layers
        for i in range(self.nlayers):
            layername = f'fc_{i+1}'
            layermodule = nn.Linear(self.layersizes[i], self.layersizes[i+1])
            self.layernames.append(layername)
            setattr(self, layername, layermodule)

    def forward(self, x):
        """
        Implement the forward pass of the NN.
        """
        for i, layername in enumerate(self.layernames):
            fclayer = getattr(self, layername)
            x = fclayer(x)
            if i == self.nlayers - 1:
                x = self.final_layer_activation(x)
            else:
                x = self.activation(x)
        return x

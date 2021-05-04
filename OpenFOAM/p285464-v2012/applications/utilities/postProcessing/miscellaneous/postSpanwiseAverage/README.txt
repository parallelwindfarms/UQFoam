#-----------------------------------------------------------------------------#
# Tool to perform averaging in one homogeneous direction of an input field
# author: Marian Fuchs   
# contact: marian.fuchs@cfd.tu-berlin.de
# version: 2012_11_15  
#-----------------------------------------------------------------------------#

info:
- tool is based on the OpenFOAM utility 'postChannel'
- to compile the utility: source code should be copied to
> WM_PROJECT_DIR/applications/utilities/postProcessing/miscellaneous/
> 'wmake' to compile the source code

how-to-run:
- define seed patch and axis directions in 'postSpanwiseAverageDict', which is placed
  in the 'constant' dictionary of the current case
- run the utility by typing:
> postSpanwiseAverage
- averaged fields are saved to:
> caseDictionary/span-averaged_field

remarks:
- the utility is only applicable to structured meshes

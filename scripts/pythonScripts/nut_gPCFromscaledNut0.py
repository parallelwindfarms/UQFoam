#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 19:02:54 2020

@author: jigar
"""


import pyvista as vtki

## grid is the central object in VTK where every field is added on to grid
grid = vtki.UnstructuredGrid('./VTK/baseline_StSt_nutPCE_10000.vtk')

## cell and point-wise information of geometry is contained
cells = grid.cells
pts   = grid.points

## get a dictionary contains all cell/point information
# note that cell-based and point-based are in different size
cells_array =  grid.cell_arrays 
pts_array   = grid.point_arrays

## get a field in numpy array
nut = cells_array['nut']

## create a new cell field of pressure^2
nut1 = nut/10
nut2 = nut/100
grid._add_cell_array(nut1, 'nut1')
grid._add_cell_array(nut2, 'nut2')

## remember to save the modified vtk
grid.save('./VTK/baseline_StSt_nutPCE_10000_nutPCE.vtk')
import numpy as np
import pyvista as pv

# %% Mesh data
def getMeshData(vtkFile, D, h, ADloc, ny, nz, d_D, x_up_D_WT):
    mesh    = pv.UnstructuredGrid(vtkFile)
    nCells  = mesh.n_cells
    cCenter = np.array(mesh.cell_centers().points)
    
    xStart_WT, xEnd_WT = ADloc[0]-x_up_D_WT*D, ADloc[0]+(d_D-x_up_D_WT)*D
    cCenterWT_idx = np.where(
        (cCenter[:,0] >= xStart_WT) & (cCenter[:,0] <= xEnd_WT)
    )[0]
    cCenter_WT = cCenter[cCenterWT_idx]
    nCells_WT = cCenter_WT.shape[0]
    nx_WT = int(nCells_WT/ny/nz)
    assert nx_WT==nz==ny, f'nx_WT,nz,ny = {nx_WT,nz,ny}'
    mlMeshShape = tuple([nz, ny, nx_WT])
    
    startPlane_WT_idx = (cCenter_WT[:,0] == cCenter_WT[:,0].min())
    endPlane_WT_idx   = (cCenter_WT[:,0] == cCenter_WT[:,0].max())
    y0Plane_WT_idx = (cCenter_WT[:,1] > np.abs(cCenter_WT[:,1]).min()-1) * \
                     (cCenter_WT[:,1] < np.abs(cCenter_WT[:,1]).min()+1)
    zhPlane_WT_idx = ((cCenter_WT[:,2] > (h)) * (cCenter_WT[:,2] < (h+3)))
    
    
    cellsInDiskAtHubHeight = np.array(np.where(
        (cCenter_WT[:,0]>=xStart_WT) & (cCenter_WT[:,0]<=(xStart_WT+10))  &
        (cCenter_WT[:,1]>=-0.5*D) & (cCenter_WT[:,1]<=0.5*D) &
        (cCenter_WT[:,2]>=(h-D/2.)) & (cCenter_WT[:,2]<=(h+D/2.)) ))[0]
    
    return mesh, nCells, mlMeshShape, nCells_WT, cCenter_WT, \
        cCenterWT_idx, startPlane_WT_idx, endPlane_WT_idx, \
        y0Plane_WT_idx, zhPlane_WT_idx, cellsInDiskAtHubHeight

# %% Basecase data
def getCaseData(myUQlib, mesh, nCells_WT, cCenterWT_idx, 
                    cellsInDiskAtHubHeight, Uref):
    UMag = np.linalg.norm(mesh.cell_data['U'][cCenterWT_idx], axis=1)
    UHub = UMag[cellsInDiskAtHubHeight].mean()
    defU = (UHub-UMag)/UHub

    tke = np.array(mesh.cell_data['k'][cCenterWT_idx])
    TI = np.sqrt(tke*2/3)/UHub *100
    TIHub = TI[cellsInDiskAtHubHeight].mean()

    R = mesh.cell_data['turbulenceProperties:R'][cCenterWT_idx]
    A = myUQlib.anisotropyTensor(
        myUQlib.symmTensorToTensorv2021(R, nCells_WT), tke, nCells_WT
    )
    A = myUQlib.getSymmTensorValues(A, nCells_WT)

    nut = np.array(mesh.cell_data['nut'][cCenterWT_idx])

    return  UMag, UHub, defU, tke, TI, TIHub, R, A, nut
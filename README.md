# MeshVoxelizer.jl
A mesh voxelization function for use with Julia. Takes a GeometryBasics mesh and converts it to a boolean array. Used for meshing finite-difference grids from STL files. 

## Dependencies
* **Suggested**
    * [MeshIO/FileIO](https://juliapackages.com/p/meshio)
* **Required**
    * [GeometryBasics](https://www.juliapackages.com/p/geometrybasics)

## Usage
```Julia
using Main.Voxelizer
using FileIO
using GLMakie
#Load using MeshIO
mesh = load(path);
#Voxelize the actual mesh
@time grid = Voxelizer.voxelizeMesh(mesh,resolution=100);
#Display as a GLMakie volume
display(volume(grid,isorange = 1, isovalue=1, algorithm=:iso))
```
Optional arguments:

* `resolution` Defines the length of the largest dimension of the voxel array. Default = 100.
* `padMesh` Pads the voxel array with 1 additional element per side, offsets the elements by 1.
 
Resolution ranging from 10-500 as presented in `examples.ipynb`

![](https://i.imgur.com/tcvFloL.gif)
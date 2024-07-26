# Fite2dmt.jl
_Finite Triangular Elements for Two-Dimensional Magnetotelluric_

## Description
Fite2dmt.jl is a Julia module for two-dimensional magnetotelluric inversion that works for unstructured triangular meshes. The unstructured mesh is handy for dealing with complex geometries, such as the undulation of surface topography that can sometimes influence MT data. The effect of topography to the MT data can be reproduced by representing the topography accurately in the inversion using an unstructured mesh. In addition, the unstructured mesh can accommodate local refinement of grids at certain regions to improve modeling accuracy, such as around observation stations.

## Features

The code implements:
- Node-based finite element formulation for triangular element
- Regularized least-squares objective function
- Data-space Gauss-Newton minimization algorithm
- Reciprocity method for assembling the sensitivity matrix

The code can invert:
- apparent resistivity
- phase
- vertical magnetic field transfer function (or Tipper)

## Example

![ForGithub_2](https://user-images.githubusercontent.com/65894100/201507763-0807b98d-54d5-4545-abb7-6a51b1a88332.png)

## Documentation

Diba D, Nurhasan, Uyeshima M, Usui Y (2024) Two-dimensional magnetotelluric inversion using unstructured triangular mesh implemented in Julia. Journal of Physics: Conference Series, **2734** 012008. https://doi.org/10.1088/1742-6596/2734/1/012008

## Terms Of Use

Please cite the original documentation of FITE2DMT when publishing any work that uses this program. Proper citation acknowledges the authors' contributions and helps support further development.

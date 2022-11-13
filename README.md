# FITE2DMT
FInite Triangular Elements for Two-Dimensional MagnetoTelluric

## Description
FITE2DMT is a two-dimensional finite element magnetotelluric inversion code that works for unstructured triangular electrical resistivity structures. The unstructured mesh is handy for dealing with complex geometries. For example, the undulation of surface topography can sometimes influence MT data. This topographic effect must be corrected by accurately modeling the topography, or it might lead to an incorrect interpretation of the structure. For that, the unstructured mesh is preferred over the simple structured mesh. Besides, the unstructured mesh can accommodate local refinement to improve modeling accuracy. Common in practice is the refinement of the elements surrounding observation stations.

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
- inter-station horizontal magnetic field transfer function

## Example

![ForGithub_2](https://user-images.githubusercontent.com/65894100/201507763-0807b98d-54d5-4545-abb7-6a51b1a88332.png)

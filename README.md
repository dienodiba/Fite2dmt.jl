# MTri
Magnetotelluric Inversion with Unstructured Triangular Mesh 

## Description
MTri is a two-dimensional finite element magnetotelluric inversion code that supports unstructured mesh of triangular electrical resistivity elements. The unstructured mesh is very useful for dealing with complex geometries. Surface topography can be represented accurately so that misinterpretation due to the topographic effect can be avoided. Also, the mesh around observation stations, or wherever necessary, can be refined locally.

## Features

The code implements:
- Node-based finite element formulation for triangular elements
- Regularized least-squares objective function
- Model space Gauss-Newton minimization algorithm
- Reciprocity method for assembling the sensitivity matrix

The code can invert:
- apparent resistivity
- phase
- vertical magnetic field transfer function
- inter-station horizontal magnetic field transfer function

## Example

![ForGithub](https://user-images.githubusercontent.com/65894100/196236640-85baaff2-f9d6-4a2f-b547-07a57d6d1342.png)

# USTRIM
UnStructured TRIangular Mesh Magnetotelluric Inversion

## Description
USTRIM is a two-dimensional finite element magnetotelluric inversion code. The resistivity is discretized using unstructured mesh of triangular elements, so the topography or other irregular structures can be represented accurately.   

## Features

This code uses:
- Node-based finite element formulation on triangular elements
- Regularized least-squares objective function
- Model space Gauss-Newton minimization algorithm
- Reciprocity method for assembling the sensitivity matrix

This code is able to invert:
- apparent resistivity
- phase
- vertical magnetic field transfer function
- inter-station horizontal magnetic field transfer function

## Example

![ForGithub](https://user-images.githubusercontent.com/65894100/196236640-85baaff2-f9d6-4a2f-b547-07a57d6d1342.png)

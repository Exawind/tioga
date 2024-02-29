# TIOGA - Topology Independent Overset Grid Assembler

[![TIOGA CI](https://github.com/Exawind/tioga/actions/workflows/ci.yml/badge.svg)](https://github.com/Exawind/tioga/actions/workflows/ci.yml)

Tioga is a library for overset grid assembly on parallel distributed systems
Copyright (C) 2015 Jay Sitaraman

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

Contact: Jay Sitaraman, jsitaraman@gmail.com, (301) 741 3216, Parallel Geometric Algorithms LLC, 814 S Mary Ave, Sunnyvale, CA 94087

## Functionality

TIOGA can perform overset grid connectivity in 3-D between
multiple unstructured (or structured) meshes that are in a distributed
computing environment, i.e. each mesh is partitioned in to multiple
parts. It can accept high-order call-back functions to perform p-consistent
interpolation and searches for formulations with internal degrees of freedom
within a computational element.

## Notes

TIOGA is free software since it was developed in the personal
time of the author. It is expected to serve as an academic/research
counterpart for PUNDIT (which is the product of the CREATE A/V program
and is export controlled). TIOGA has a subset of the functionality of
PUNDIT and is about 3x slower owing to the use of Alternating Digital Tree (ADT)
searches as the baseline point-location algorithm.

## News

TIOGA is currently under development with funding from the DoE ExaWind program
towards developing overset capability in the DoE ExaWind codes for large scale
wind farm simulations.

## References

Michael J. Brazell, Jayanarayanan Sitaraman, Dimitri J. Mavriplis, An overset mesh approach for 3D mixed element
high-order discretizations, In Journal of Computational Physics, Volume 322, 2016,
Pages 33-51, ISSN 0021-9991, https://doi.org/10.1016/j.jcp.2016.06.031.
(http://www.sciencedirect.com/science/article/pii/S002199911630256X)


Roget, B. and Sitaraman, J., "Robust and efficient overset grid assembly for partitioned unstructured meshes",
Journal of Computational Physics, v 260, March 2014, Pages 1-24

Brazell, M., Sitaraman J. and Mavriplis D.,"An Overset Mesh Approach for 3D Mixed
Element High Order Discretizations", Proceedings of 2014 Overset Grid Symposium,
Atlanta, GA, Oct 6-9, 2014.
http://2014.oversetgridsymposium.org/assets/presentations/3_1/Brazell_ogs_2014.pdf

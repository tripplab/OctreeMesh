// Shape.h
// Copyright (C) 2010 Miguel Vargas (miguel.vargas@gmail.com)
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef _Shape_h_
#define _Shape_h_


enum ShapeType
{
	shape_undefined     = 0,
	shape_linear        = 1,
	shape_triangle      = 2,
	shape_quadrilateral = 3,
	shape_tetrahedron   = 4,
	shape_hexahedron    = 5
};


#define FaceTypeMacro(element_type) (((element_type) == shape_triangle || (element_type) == shape_quadrilateral) ? shape_linear : ((element_type) == shape_tetrahedron) ? shape_triangle : ((element_type) == shape_hexahedron) ? shape_quadrilateral : shape_undefined)
#define NodesPerFacetMacro(element_type, nodes_per_element) ((element_type) == shape_triangle ? ((nodes_per_element) == 3 ? 2 : 3) : (element_type) == shape_quadrilateral ? ((nodes_per_element) == 4 ? 2 : 3) : (element_type) == shape_tetrahedron ? ((nodes_per_element) == 4 ? 3 : 6) : (element_type) == shape_hexahedron ? ((nodes_per_element) == 8 ? 4 : ((nodes_per_element) == 20) ? 8 : 9) : 0)
#define FacetsPerElement(element_type) ((element_type) == shape_triangle ? 3 : (element_type) == shape_quadrilateral ? 4 : (element_type) == shape_tetrahedron ? 4 : (element_type) == shape_hexahedron ? 6 : 0)

typedef int Facet[4];
extern const Facet* facets[];

#define FacetNode(element_type, face, node) (facets[(element_type) - 1][(face) - 1][(node) - 1])

#endif

// Shape.cpp
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

#include <FiniteElement/Shape.h>


const Facet linear_facets[1] = {{1, 2}};
const Facet triangle_facets[3] = {{1, 2}, {2, 3}, {3, 1}};
const Facet quadrilateral_facets[4] = {{1, 2}, {2, 3}, {3, 4}, {4, 1}};
const Facet tetrahedra_facets[4] = {{1, 2, 3}, {2, 4, 3}, {3, 4, 1}, {4, 2, 1}};
const Facet hexahedra_facets[6] = {{1, 2, 3, 4}, {1, 4, 8, 5}, {1, 5, 6, 2}, {2, 6, 7, 3}, {3, 7, 8, 4}, {5, 8, 7, 6}};
const Facet* facets[] = {linear_facets, triangle_facets, quadrilateral_facets, tetrahedra_facets, hexahedra_facets};

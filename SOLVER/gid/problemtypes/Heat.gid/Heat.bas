*intformat "%i "
*realformat "%.7g "
*if((nelem(Point)>0)||(nelem(Linear)>0)||(nelem(Prism)>0)||(nelem(Pyramid)>0)||(nelem(Sphere)>0))
*messagebox Error: Element types supported: Triangle, Quadrilateral, Tetrahedra, Hexahedra
*endif
*if(ndime==2)
*if((nelem(Triangle)>0)&&(nelem(Quadrilateral)>0))
*messagebox Error: Only single element type meshes are supported.
*endif
*if(nelem(Triangle)>0)
*set var elementtype=2
*else
*set var elementtype=3
*endif
*endif
*if(ndime==3)
*if((nelem(Tetrahedra)>0)&&(nelem(Hexahedra)>0))
*messagebox Error: Only single element type meshes are supported.
*endif
*if(nelem(Tetrahedra)>0)
*set var elementtype=4
*else
*set var elementtype=5
*endif
*endif
; Geometry file

{Nodes}
*ndime; Dimension
*set elems(all)
*npoin; Nodes count
; X1 X2 ...
*loop nodes
*nodescoord
*end nodes

{Mesh}
*elementtype; Element type (2=Triangle, 3=Quadrilateral, 4=Tetrahedra, 5=Hexahedra)
*nnode; Nodes per element
*nelem; Elements count
; Material Node1 Node2 ...
*loop elems
*elemsmat*elemsconec
*end elems

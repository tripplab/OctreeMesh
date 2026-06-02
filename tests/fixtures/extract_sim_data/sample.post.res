GiD Post Result File 1.0

GaussPoints "StrainStressGaussPoints" ElemType Hexahedra
Number Of Gauss Points: 8
Natural Coordinates: Given
 -5.77350e-01 -5.77350e-01 -5.77350e-01
 -5.77350e-01 -5.77350e-01 5.77350e-01
 -5.77350e-01 5.77350e-01 -5.77350e-01
 -5.77350e-01 5.77350e-01 5.77350e-01
 5.77350e-01 -5.77350e-01 -5.77350e-01
 5.77350e-01 -5.77350e-01 5.77350e-01
 5.77350e-01 5.77350e-01 -5.77350e-01
 5.77350e-01 5.77350e-01 5.77350e-01
End GaussPoints

Result "Displacement" "Solid" 1 Vector OnNodes
Values
1 1.0 0.0 0.0
2 0.0 2.0 0.0
3 0.0 0.0 3.0
End Values

Result "Von Mises" "Solid" 1 Scalar OnGaussPoints "StrainStressGaussPoints"
Values
1 1.0
2.0
3.0
4.0
5.0
6.0
7.0
8.0
2 9.0 10.0
11.0 12.0
13.0
14.0
15.0
16.0
End Values

*intformat "%i "
*realformat "%.7g "
*if(strcmp(gendata(Plane_problem),"PlaneStress")==0)
*set var planeproblem=1
*elseif(strcmp(gendata(Plane_problem),"PlaneStrain")==0)
*set var planeproblem=2
*endif
*if(strcmp(gendata(Problem_type),"Stationary")==0)
*set var problemtype=1
*elseif(strcmp(gendata(Problem_type),"Dynamic")==0)
*set var problemtype=2
*endif
; Problem file

{General}
*planeproblem; Plane_problem (1=PlaneStress, 2=PlaneStrain)
*gendata(Calculate_displacement,int); Calculate_displacement (0=No, 1=Yes)
*gendata(Calculate_strain,int); Calculate_strain (0=No, 1=Yes)
*gendata(Calculate_stress,int); Calculate_stress (0=No, 1=Yes)
*gendata(Calculate_Von_Mises,int); Calculate_Von_Mises (0=No, 1=Yes)
*problemtype; Problem_type (1=Stationary, 2=Dynamic)
*gendata(Save_mesh,int); Save_mesh (0=No, 1=Yes)
*gendata(Save_system_of_equations,int); Save_system_of_equations (0=No, 1=Yes)

{DynamicParameters}
*gendata(Time_per_step,real); Time_per_step
*gendata(Steps,int); Steps
*gendata(Result_every_steps,int); Result_every_steps
*gendata(Time_scheme_factor,real); Time_scheme_factor (The alpha-mehtod is second-order accurate and unconditionally stable for [-1/3, 0])
*gendata(Rayleigh_damping_a,real); Rayleigh damping coefficient a
*gendata(Rayleigh_damping_b,real); Rayleigh damping coefficient b

{Gravity}
*gendata(Use_mass_forces,int); Use_mass_forces (0=No, 1=Yes)
*gendata(Gravity,real); Gravity

{UserFunctions}
; f1, f2, ..., f20 (one per line)
*gendata(f1)
*gendata(f2)
*gendata(f3)
*gendata(f4)
*gendata(f5)
*gendata(f6)
*gendata(f7)
*gendata(f8)
*gendata(f9)
*gendata(f10)
*gendata(f11)
*gendata(f12)
*gendata(f13)
*gendata(f14)
*gendata(f15)
*gendata(f16)
*gendata(f17)
*gendata(f18)
*gendata(f19)
*gendata(f20)

{Materials}
*nmats; Materials count
; Poisson_ratio Young_modulus Density Thickness
*loop materials
*matprop(Poisson_ratio,real)*matprop(Young_modulus,real)*matprop(Density,real)*matprop(Thickness,real); *matprop(0)
*end materials

{Displacement}
*set cond Displacement_Point *nodes *canrepeat
*add cond Displacement_Line *nodes
*add cond Displacement_Surface *nodes
*condnumentities(int); Displacements count
*if(*condnumentities(int)>0)
; Node User_function_for_x Displacement_x Fixed_x User_function_for_y Displacement_y Fixed_y User_function_for_z Displacement_z Fixed_z
*loop nodes *onlyincond
*if((cond(Fixed_x,int)==0)&&(cond(Fixed_y,int)==0)&&(cond(Fixed_z,int)==0))
*nodesnum*cond(User_function_for_x,int)*cond(Displacement_x,real)*cond(Fixed_x,int)*cond(User_function_for_y,int)*cond(Displacement_y,real)*cond(Fixed_y,int)*cond(User_function_for_z,int)*cond(Displacement_z,real)*cond(Fixed_z,int)
*endif
*end nodes
*loop nodes *onlyincond
*if((cond(Fixed_x,int)!=0)||(cond(Fixed_y,int)!=0)||(cond(Fixed_z,int)!=0))
*nodesnum*cond(User_function_for_x,int)*cond(Displacement_x,real)*cond(Fixed_x,int)*cond(User_function_for_y,int)*cond(Displacement_y,real)*cond(Fixed_y,int)*cond(User_function_for_z,int)*cond(Displacement_z,real)*cond(Fixed_z,int)
*endif
*end nodes
*endif

{NodalForce}
*set cond Nodal_Force_Point *nodes
*add cond Nodal_Force_Line *nodes
*add cond Nodal_Force_Surface *nodes
*condnumentities(int); Nodal forces count
*if(*condnumentities(int)>0)
; Node User_function_for_x Force_x User_function_for_y Force_y User_function_for_z Force_z
*loop nodes *onlyincond
*nodesnum*cond(User_function_for_x,int)*cond(Force_x,real)*cond(User_function_for_y,int)*cond(Force_y,real)*cond(User_function_for_z,int)*cond(Force_z,real)
*end nodes
*endif

{NormalForce}
*set cond Normal_Force_Line *elems
*add cond Normal_Force_Surface *elems
*condnumentities(int); Normal forces count
*if(*condnumentities(int)>0)
; Element User_function Force FaceNode1 FaceNode2 ...
*loop elems *onlyincond
*elemsnum*cond(User_function,int)*cond(Force,real)*GlobalNodes
*end elems
*endif

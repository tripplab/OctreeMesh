*intformat "%i "
*realformat "%.7g "
*if(strcmp(gendata(Problem_type),"TimeIndependent")==0)
*set var problemtype=1
*elseif(strcmp(gendata(Problem_type),"TimeDependent")==0)
*set var problemtype=2
*endif
; Problem file

{General}
*gendata(Calculate_temperature,int); Calculate_temperature (0=No, 1=Yes)
*gendata(Calculate_flux,int); Calculate_flux (0=No, 1=Yes)
*problemtype; Problem_type (1=TimeIndependent, 2=TimeDependent)
*gendata(Save_mesh,int); Save_mesh (0=No, 1=Yes)
*gendata(Save_system_of_equations,int); Save_system_of_equations (0=No, 1=Yes)

{TimeDependent}
*gendata(Time_per_step,real); Time_per_step
*gendata(Steps,int); Steps
*gendata(Result_every_steps,int); Result_every_steps
*gendata(Time_scheme_factor,real); Time scheme factor (0.0=Fully explicit scheme, 0.5=Crank-Nicolson scheme, 1.0=Fully implicit scheme)

{UserFunctions}
; f1, f2, ..., f10 (one per line)
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

{Materials}
*nmats; Materials count
; Thermal_conductivity Mass_density Specific_heat_capacity
*loop materials
*matprop(Thermal_conductivity,real)*matprop(Mass_density,real)*matprop(Specific_heat_capacity,real); *matprop(0)
*end materials

{Temperature}
*set cond Temperature_Point *nodes *canrepeat
*add cond Temperature_Line *nodes
*add cond Temperature_Surface *nodes
*add cond Temperature_Volume *nodes
*condnumentities(int); Temperatures count
*if(*condnumentities(int)>0)
; Node User_function Temperature Fixed
*loop nodes *onlyincond
*if(cond(Fixed,int)==0)
*nodesnum*cond(User_function,int)*cond(Temperature,real)*cond(Fixed,int)
*endif
*end nodes
*loop nodes *onlyincond
*if(cond(Fixed,int)==1)
*nodesnum*cond(User_function,int)*cond(Temperature,real)*cond(Fixed,int)
*endif
*end nodes
*endif

{HeatFlow}
*set cond Heat_Flow_Line *elems
*add cond Heat_Flow_Surface *elems
*condnumentities(int); Heat flows count
*if(*condnumentities(int)>0)
; User_function Flux FaceNode1 FaceNode2 ...
*loop elems *onlyincond
*cond(User_function,int)*cond(Flux,real)*GlobalNodes
*end elems
*endif

{Source}
*set cond Source_Surface *elems
*add cond Source_Volume *elems
*condnumentities(int); Sources count
*if(*condnumentities(int)>0)
; Element User_function Heat
*loop elems *onlyincond
*elemsnum*cond(User_function,int)*cond(Heat,real)
*end elems
*endif

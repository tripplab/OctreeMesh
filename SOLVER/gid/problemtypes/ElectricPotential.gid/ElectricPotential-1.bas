*intformat "%i "
*realformat "%.7g "
*if(strcmp(gendata(Problem_type),"Simple")==0)
*set var problemtype=1
*elseif(strcmp(gendata(Problem_type),"Capacitance_matrix")==0)
*set var problemtype=2
*elseif(strcmp(gendata(Problem_type),"Sensitivity_analysis")==0)
*set var problemtype=3
*endif
; Problem file

{General}
*problemtype; Problem_type (1=Simple, 2=Capacitance_matrix, 3=Sensitivity_analysis)
*gendata(Save_mesh,int); Save_mesh (0=No, 1=Yes)

{Simple}
*gendata(Calculate_potential,int); Calculate_potential (0=No, 1=Yes)
*gendata(Calculate_electric_field,int); Calculate_electric_field (0=No, 1=Yes)
*gendata(Save_system_of_equations,int); Save_system_of_equations (0=No, 1=Yes)

{CapacitanceMatrix}
*gendata(Electrodes_count,int); Electrodes_count
*gendata(Electrodes_voltage,real); Electrodes_voltage

{SensitivityAnalysis}
*gendata(Electrode_segments,int); Electrode_segments
*gendata(Segments_per_step,int); Segments_per_step
*gendata(Result_on_nodes,int); Result_on_nodes (0=No, 1=Yes)

{Materials}
*operation(nmats+2); Materials count
; Sensitivity_use Relative_permittivity
*loop materials
*matprop(Sensitivity_use,int)*matprop(Relative_permittivity,real); *matprop(0)
*end materials
0 *gendata(High_permittivity,real); Sensitivity_analysis_high_permittivity
0 *gendata(Low_permittivity,real); Sensitivity_analysis_low_permittivity

{Potential}
*set cond Electric_potential_Point *nodes
*add cond Electric_potential_Line *nodes
*add cond Electric_potential_Surface *nodes
*add cond Electric_potential_Volume *nodes
*condnumentities(int); Potentials count
*if(*condnumentities(int)>0)
; Node Potential Electrode
*loop nodes *onlyincond
*nodesnum*cond(Potential,real)*cond(Electrode,int)
*end nodes
*endif

{ElectricField}
*set cond Electric_field_Line *elems
*add cond Electric_field_Surface *elems
*condnumentities(int); Electric_field count
*if(*condnumentities(int)>0)
; Field_strength FaceNode1 FaceNode2 ...
*loop elems *onlyincond
*cond(Field_strength,real)*GlobalNodes
*end elems
*endif

{Source}
*set cond Charge_Surface *elems
*add cond Charge_Volume *elems
*condnumentities(int); Sources count
*if(*condnumentities(int)>0)
; Element Charge_density
*loop elems *onlyincond
*elemsnum*cond(Charge_density,real)
*end elems
*endif

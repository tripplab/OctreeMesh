#!/usr/bin/awk -f

## by trippm@bmd [tripplab.com] 02/02/2022

## Bash script to rotate a set of atom positions in PDB format according to a rotation matrix
## ./apply-matrix.awk original.pdb rotate.mtx > rotated.pdb

BEGIN{

#COLUMNS        DATA TYPE       FIELD         DEFINITION
#---------------------------------------------------------------------------
# 1 -  6        Record name     "ATOM  "
# 7 - 11        Integer         serial        Atom serial number.
#13 - 16        Atom            name          Atom name.
#17             Character       altLoc        Alternate location indicator.
#18 - 20        Residue name    resName       Residue name.
#22             Character       chainID       Chain identifier.
#23 - 26        Integer         resSeq        Residue sequence number.
#27             AChar           iCode         Code for insertion of residues.
#31 - 38        Real(8.3)       x             Orthogonal coordinates for X in
#                                             Angstroms.
#39 - 46        Real(8.3)       y             Orthogonal coordinates for Y in
#                                             Angstroms.
#47 - 54        Real(8.3)       z             Orthogonal coordinates for Z in
#                                             Angstroms.
#55 - 60        Real(6.2)       occupancy     Occupancy.
#61 - 66        Real(6.2)       tempFactor    Temperature factor.
#73 - 76        LString(4)      segID         Segment identifier,left-justified
#77 - 78        LString(2)      element       Element symbol,right-justified
#79 - 80        LString(2)      charge        Charge on the atom.

# topp format
#  getline<ARGV[2];
  for(i=1;i<=3;i++) {
    getline<ARGV[2];
    rot[i,1]=$1;
    rot[i,2]=$2;
    rot[i,3]=$3;
#    print "r",rot[i,1],rot[i,2],rot[i,3];

# stamp format
#    trans[i]=$4;
  }

# topp format
  getline<ARGV[2];
  trans[1]=0; #$1;
  trans[2]=0; #$2;
  trans[3]=0; #$3;
#  print "t",trans[1],trans[2],trans[3];

  lines=0;
  while((getline<ARGV[1])>0){

   ##if($1=="ATOM") {
   if( ($1~"ATOM") && (length($4)==3) && ($4!="SOL") ) {
    lines++;

    rec_name="";
    serial="";
    name="";
    altLoc="";
    resName="";
    chainID="";
    resSeq="";
    iCode="";
    n_char=split($0,char,"");
    for(i=1;i<=6;i++) rec_name=rec_name char[i];
    for(i=7;i<=11;i++) serial=serial char[i];
    for(i=13;i<=16;i++) name=name char[i];
    name=$3;
    altLoc= char[17];
    for(i=18;i<=20;i++) resName=resName char[i];
    chainID= char[22];
    for(i=23;i<=26;i++) resSeq=resSeq char[i];
    iCode= char[27];

#    rec_init="";
    x="";
    y="";
    z="";
    occupancy="";
    tempFactor="";
    segID="";
    element="";
    charge=""
#    n_char=split($0,char,"");
#    for(i=1;i<=30;i++) rec_init=rec_init char[i];
    for(i=31;i<=38;i++) x=x char[i];
    for(i=39;i<=46;i++) y=y char[i];
    for(i=47;i<=54;i++) z=z char[i];

    xrot=(rot[1,1]*x)+(rot[1,2]*y)+(rot[1,3]*z)
    yrot=(rot[2,1]*x)+(rot[2,2]*y)+(rot[2,3]*z)
    zrot=(rot[3,1]*x)+(rot[3,2]*y)+(rot[3,3]*z)
    x=trans[1]+xrot
    y=trans[2]+yrot
    z=trans[3]+zrot

    for(i=55;i<=60;i++) occupancy=occupancy char[i];
    for(i=61;i<=66;i++) tempFactor=tempFactor char[i];
    for(i=73;i<=76;i++) segID=segID char[i];
    for(i=77;i<=78;i++) element=element char[i];
    for(i=79;i<=80;i++) charge=charge char[i];

#    printf("%30s%8.3f%8.3f%8.3f%6.2f%6.2f       %4s%2s%2s\n",rec_init,x,y,z,occupancy,tempFactor,segID,element,charge);
    printf("%6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f       %4s%-2s%2s\n",rec_name,serial,name,altLoc,resName,chainID,resSeq,iCode,x,y,z,occupancy,tempFactor,segID,element,charge);
   }
  }
#print "Trasnformed ",lines," lines using matrix in ",ARGV[2];
}

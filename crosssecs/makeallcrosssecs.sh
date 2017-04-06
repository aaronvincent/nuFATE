#!/bin/bash

#  makeallcrosssecs.sh
#  This just runs the resulting code from the compiled fortran files.
#
#  Created by Aaron Vincent on 2017-03-06.
#
Emin=3
Emax=10
#Total cross sections
./CCxstot.x $Emin $Emax -1 ../../data/nuebarxs.dat
./CCxstot.x $Emin $Emax -2 ../../data/numubarxs.dat
./CCxstot.x $Emin $Emax -3 ../../data/nutaubarxs.dat
./CCxstot.x $Emin $Emax 1 ../../data/nuexs.dat
./CCxstot.x $Emin $Emax 2 ../../data/numuxs.dat
./CCxstot.x $Emin $Emax 3 ../../data/nutauxs.dat

#Neutral current: the same for all flavours
./CCxsdiff.x $Emin $Emax -1 ../../data/dxsnubar.dat
./CCxsdiff.x $Emin $Emax 1 ../../data/dxsnu.dat

#Tau regeneration.
./tauRegen.x $Emin $Emax -3 ../../data/tbarfull.dat
./tauRegen.x $Emin $Emax 3 ../../data/tfull.dat

#Secondary production. e and mu are the same.
#./tauRegen.x $Emin $Emax -2 ../../data/secbarfull.dat
#./tauRegen.x $Emin $Emax 2 ../../data/secfull.dat


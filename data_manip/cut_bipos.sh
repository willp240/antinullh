source /home/kroupova/env_git_rat.sh
cd /home/kroupova/bb_sigex/data_manip

g++ remove_notrigs.cpp $(root-config --cflags --libs) -o remove_notrigs -Iinclude -I/home/kroupova/rat/include -Iinclude/RAT -I/home/kroupova/rat/include/RAT -I/home/kroupova/rat/include/libpq -I/home/kroupova/rat/include/RAT/DS -I/data/snoplus/software/snocave_SL6/clhep-2.1.1.0/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include/Geant4 -L/home/kroupova/rat/lib -lRATEvent_Linux -lrat_Linux
g++ cut_bipos.cpp $(root-config --cflags --libs) -o cut_bipos -Iinclude -I/home/kroupova/rat/include -Iinclude/RAT -I/home/kroupova/rat/include/RAT -I/home/kroupova/rat/include/libpq -I/home/kroupova/rat/include/RAT/DS -I/data/snoplus/software/snocave_SL6/clhep-2.1.1.0/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include/Geant4 -L/home/kroupova/rat/lib -lRATEvent_Linux -lrat_Linux

TESTING=0
DELTA_T1=500
DELTA_T2=3936000
DELTA_R=1500

./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo212/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo212_triggeredonly/trigs_only.root $TESTING

./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo214/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo214_triggeredonly/trigs_only.root $TESTING



./cut_bipos /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo212_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo212_cleaned/cleaned.root $DELTA_T1 $DELTA_T2 $DELTA_R /home/kroupova/bb_sigex/data_manip/bipo_cut_eff/Bipo212 $TESTING

./cut_bipos /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo214_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedBipo214_cleaned/cleaned.root $DELTA_T1 $DELTA_T2 $DELTA_R /home/kroupova/bb_sigex/data_manip/bipo_cut_eff/Bipo214 $TESTING




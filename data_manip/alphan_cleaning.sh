source /home/kroupova/env_git_rat.sh
cd /home/kroupova/bb_sigex/data_manip

g++ remove_notrigs.cpp $(root-config --cflags --libs) -o remove_notrigs -Iinclude -I/home/kroupova/rat/include -Iinclude/RAT -I/home/kroupova/rat/include/RAT -I/home/kroupova/rat/include/libpq -I/home/kroupova/rat/include/RAT/DS -I/data/snoplus/software/snocave_SL6/clhep-2.1.1.0/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include/Geant4 -L/home/kroupova/rat/lib -lRATEvent_Linux -lrat_Linux
g++ clean_alpha_n.cpp $(root-config --cflags --libs) -o clean_alpha_n -Iinclude -I/home/kroupova/rat/include -Iinclude/RAT -I/home/kroupova/rat/include/RAT -I/home/kroupova/rat/include/libpq -I/home/kroupova/rat/include/RAT/DS -I/data/snoplus/software/snocave_SL6/clhep-2.1.1.0/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include -I/data/snoplus/software/snocave_SL6/geant4.10.0.p02/include/Geant4 -L/home/kroupova/rat/lib -lRATEvent_Linux -lrat_Linux

TESTING=0
EMIN=1
DELTA_T=1000000

./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_13c/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_13c_triggeredonly/trigs_only.root $TESTING
./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Av_13c/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Av_13c_triggeredonly/trigs_only.root $TESTING
./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Ls_13c/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Ls_13c_triggeredonly/trigs_only.root $TESTING
./remove_notrigs "/data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avout_Av_13c/*.root" /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avout_Av_13c_triggeredonly/trigs_only.root $TESTING


./clean_alpha_n /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_13c_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_13c_cleaned/cleaned.root $EMIN $DELTA_T /home/kroupova/bb_sigex/data_manip/alphan_cut_eff/Alphan_Telab_13c $TESTING

./clean_alpha_n /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Av_13c_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Av_13c_cleaned/cleaned.root $EMIN $DELTA_T /home/kroupova/bb_sigex/data_manip/alphan_cut_eff/Alphan_Telab_Avin_Av_13c $TESTING

./clean_alpha_n /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Ls_13c_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avin_Ls_13c_cleaned/cleaned.root $EMIN $DELTA_T /home/kroupova/bb_sigex/data_manip/alphan_cut_eff/Alphan_Telab_Avin_Ls_13c $TESTING

./clean_alpha_n /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avout_Av_13c_triggeredonly/trigs_only.root /data/snoplus/griddata/Prod_Rat6163_TeLoaded/TeLoadedAlphan_Telab_Avout_Av_13c_cleaned/cleaned.root $EMIN $DELTA_T /home/kroupova/bb_sigex/data_manip/alphan_cut_eff/Alphan_Telab_Avout_Av_13c $TESTING

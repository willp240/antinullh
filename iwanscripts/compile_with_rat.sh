g++ $1 -o $(basename "$1" .cpp) -O2 -I/data/snoplus/software/snocave_CentOS7/root-5.34.36/include -I/home/lidgard/rat-7/include -Iinclude/RAT -I/home/lidgard/rat-7/include/RAT -I/home/lidgard/rat-7/include/libpq -I/home/lidgard/rat-7/include/RAT/DS -I/home/lidgard/rat-7/include/RAT/DU -I/data/snoplus/software/snocave_CentOS7/gsl-1.16/include -I/home/lidgard/oxsx/include -L/data/snoplus/software/snocave_CentOS7/root-5.34.36/lib -I/data/snoplus/software/snocave_CentOS7/clhep/x86_64-cc7-gcc48-opt/include -I/data/snoplus/software/snocave_CentOS7/geant4.10.0.p02/include -L/data/snoplus/software/snocave_CentOS7/geant4.10.0.p02/lib64 -L/data/snoplus/software/snocave_CentOS7/gsl-1.16/lib -L/home/lidgard/oxsx/build -loxsx -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -lG4run -lG4clhep -lCore -lRIO -lHist -lGpad -lTree -lRint -lPostscript -lMatrix -lMathCore -lThread -lpthread -lMinuit2 -lpthread -lm -ldl -lMinuit2 -lGraf -larmadillo -lgsl -lgslcblas -lgsl -lgslcblas -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -larmadillo -Wl,-rpath,/data/snoplus/software/snocave_CentOS7/root-5.34.36/lib,-rpath,/data/snoplus/software/snocave_CentOS7/gsl-1.16/lib -L/home/lidgard/rat-7/lib -lRATEvent_Linux -lrat_Linux
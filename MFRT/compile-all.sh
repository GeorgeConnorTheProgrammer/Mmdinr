#!/bin/bash

# Use GNU to get mfrt.o
g++ -Wall -I/$HOME/apps/root/include -c mfrt.cpp

# Use mfrt.o to link libraries and compile using GNU into exec format
g++ -O2 -m64 mfrt.o -L/$HOME/apps/root/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -o mfrtROOT -L/$HOME/apps/root/lib -Wl,-R/$HOME/apps/root/lib

# Clean up
rm -f ./mfrt.o

exit 0

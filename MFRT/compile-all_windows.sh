#!/bin/bash

# Use GNU to get mfrt.o
g++ -Wall -I<path to root include, example: /../root/include> -c mfrt.cpp

# Use mfrt.o to link libraries and compile into exec format
g++ -O2 -m64 mfrt.o -L</../root/lib> -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -o mfrtROOT -L</../root/lib> -Wl,-R</../root/lib>

# Clean up
rm -f ./mfrt.o

exit 0

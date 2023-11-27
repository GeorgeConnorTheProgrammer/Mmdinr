#!/bin/bash

g++ mfrt.cpp $(root-config --glibs --cflags --libs) -o mfrt

exit 0

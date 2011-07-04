#!/bin/bash
rm RunTest
rm *.o
rm *.hi
ghc --make RunTest.hs -rtsopts -threaded -O2 -o RunTest


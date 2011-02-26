#!/bin/bash
rm RunTest
rm *.o
rm *.hi
ghc --make RunTest.hs -O2 -o RunTest


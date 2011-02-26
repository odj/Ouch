rm *.hi *.o
rm mapcount
ghc --make mapcount.hs -O2 -o mapcount


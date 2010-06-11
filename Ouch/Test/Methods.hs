{-------------------------------------------------------------------------------
-------------------------------------------------------------------------------
    Methods - a module to define unit tests
    
    Copyright (c) 2010 Orion D. Jankowski
    
    This file is part of Ouch, a chemical informatics toolkit
    written entirely in Haskell.

    Ouch is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Ouch is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Ouch.  If not, see <http://www.gnu.org/licenses/>.

-------------------------------------------------------------------------------
------------------------------------------------------------------------------}

module Ouch.Test.Methods 
    (
       TestData(..)
     , tests
     , performTests
     , testTest
     , testFail
     , testMolForm
     ) where
         
{-# LANGUAGE RecordWildCards, CPP #-}
import Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Input.Smiles
import Ouch.Data.Atom

import Data.Either
import Data.Maybe


-- Data structure to hold test/result pairs and descriptions.  All 'functions' should
-- evaluate as (String -> Either String String) 
data TestData = TestData { function :: (String -> Either String String)  
                         , description :: String  
                         , input :: String
                         , outcome :: String    
                         } deriving (Show)

performTests :: [TestData] -> String
performTests [] = ""
performTests (x:xs) = output
 where result = (function x) (input x)
       output = case result of
           Left s ->  "\n" ++ description x ++ ":\t" ++ "FAILED\\t-with error string:\t" ++ s ++ "\n" ++ performTests xs
           Right s -> if s == (outcome x)
                      then "." ++ performTests xs
                      else "\n" ++ description x ++ ":\t" ++ "FAILED\t-with output: " ++ s ++ " || should get: " ++ (outcome x) ++ "\n" ++ performTests xs





-- Simple test function
testTest::String -> Either String String
testTest s = Right s

-- Simple test fail function 
testFail::String -> Either String String
testFail s = Left s

-- Test smiles to formula
testMolForm::String -> Either String String
testMolForm s = molecularFormula $ makeMoleculeFromSmiles s

testMolWt::String -> Either String String
testMolWt s = output
 where eitherMolWt =  molecularWeight $ makeMoleculeFromSmiles s
       output = case eitherMolWt of
           Left mw   -> Left mw
           Right mw  -> Right (show $ floor (10 * mw))



tests = [ 
           TestData {function=testTest, input="OK", description="Test SHOULD FAIL", outcome="Not OK"}
         , TestData {function=testFail, input="OK Error", description="Test SHOULD FAIL", outcome="OK"}
         , TestData {function=testTest, input="OK", description="Testing should pass", outcome="OK"}
         , TestData {function=testMolForm, input="C", description="Test Smiles-MF 1", outcome="CH4"}
         , TestData {function=testMolForm, input="CC", description="Test Smiles-MF 2", outcome="C2H6"}
         , TestData {function=testMolForm, input="CCC", description="Test Smiles-MF 3", outcome="C3H8"}
         , TestData {function=testMolForm, input="C(C)C", description="Test Smiles-MF 4", outcome="C3H8"}
         , TestData {function=testMolForm, input="C(C(C))", description="Test Smiles-MF 5", outcome="C3H8"}
         , TestData {function=testMolForm, input="C(C(C))C", description="Test Smiles-MF 6", outcome="C4H10"}
         , TestData {function=testMolForm, input="C(C(C))(C)", description="Test Smiles-MF 7", outcome="C4H10"}
         , TestData {function=testMolForm, input="C(C(C))CC(C)", description="Test Smiles-MF 8", outcome="C6H14"}
         , TestData {function=testMolForm, input="C(C(C(C)(C)))", description="Test Smiles-MF 9", outcome="C5H12"}
         , TestData {function=testMolForm, input="C(C)(C)(C)C", description="Test Smiles-MF 10", outcome="C5H12"}
         , TestData {function=testMolForm, input="C=C", description="Test Smiles-MF 11", outcome="C2H4"}
         , TestData {function=testMolForm, input="C=C=C", description="Test Smiles-MF 12", outcome="C3H4"}
         , TestData {function=testMolForm, input="C#C", description="Test Smiles-MF 13", outcome="C2H2"}
         , TestData {function=testMolForm, input="C#CC#C", description="Test Smiles-MF 14", outcome="C4H2"}
         , TestData {function=testMolForm, input="C.C", description="Test Smiles-MF 15", outcome="C2H8"}
         , TestData {function=testMolForm, input="C-C", description="Test Smiles-MF 16", outcome="C2H6"}
         , TestData {function=testMolForm, input="C-C-C-C", description="Test Smiles-MF 17", outcome="C4H10"}
         , TestData {function=testMolForm, input="C.C.C.C", description="Test Smiles-MF 18", outcome="C4H16"}
         , TestData {function=testMolForm, input="N", description="Test Smiles-MF 19", outcome="H3N"}
         , TestData {function=testMolForm, input="CN", description="Test Smiles-MF 20", outcome="CH5N"}
         , TestData {function=testMolForm, input="C#N", description="Test Smiles-MF 21", outcome="CHN"}
         , TestData {function=testMolForm, input="OCCN", description="Test Smiles-MF 22", outcome="C2H7NO"}
         , TestData {function=testMolForm, input="CBr", description="Test Smiles-MF 23", outcome="CH3Br"}
         , TestData {function=testMolForm, input="CCl", description="Test Smiles-MF 24", outcome="CH3Cl"}
         , TestData {function=testMolForm, input="CF", description="Test Smiles-MF 25", outcome="CH3F"}
         , TestData {function=testMolForm, input="OS(=O)(=O)O", description="Test Smiles-MF 26", outcome="H2O4S"}    
         , TestData {function=testMolForm, input="CC(=O)CC(=O)OC", description="Test Smiles-MF 27", outcome="C5H8O3"}
         , TestData {function=testMolForm, input="C1CC1", description="Test Smiles-MF 28", outcome="C3H6"}
         , TestData {function=testMolForm, input="C1CCCC1", description="Test Smiles-MF 29", outcome="C5H10"}
         , TestData {function=testMolForm, input="C2CCCC2", description="Test Smiles-MF 30", outcome="C5H10"}
         , TestData {function=testMolForm, input="C3CCC3", description="Test Smiles-MF 31", outcome="C4H8"}
         , TestData {function=testMolForm, input="C1CC=CC1", description="Test Smiles-MF 32", outcome="C5H8"}
         , TestData {function=testMolForm, input="C1CC1CC", description="Test Smiles-MF 33", outcome="C5H10"} 
         , TestData {function=testMolForm, input="C1.C1", description="Test Smiles-MF 34", outcome="C2H6"}
         , TestData {function=testMolForm, input="C12CCCC1CCCC2", description="Test Smiles-MF 35", outcome="C9H16"}
         , TestData {function=testMolForm, input="C123CCCC(CC5CC2)(CCC3)(CC5CC1)", description="Test Smiles-MF 36", outcome="C16H26"}   
         , TestData {function=testMolForm, input="C=1CCCCCC1", description="Test Smiles-MF 37", outcome="C7H12"}
         , TestData {function=testMolForm, input="c1ccccc1", description="Test Smiles-MF 38", outcome="C6H6"}
         , TestData {function=testMolForm, input="n1cccc1", description="Test Smiles-MF 39", outcome="C4H4N"}
         , TestData {function=testMolForm, input="cc", description="Test Smiles-MF 40", outcome="C2H4"}
         , TestData {function=testMolForm, input="o1cccc1", description="Test Smiles-MF 41", outcome="C4H4O"}
         , TestData {function=testMolForm, input="CCC(C)C(NC(=O)C(C)NC(=O)C(CC(O)=O)NC(=O)C(C)NC(=O)C(N)CC1=CC=C(O)C=C1)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(C(C)O)C(=O)NC(CC(N)=O)C(=O)NC(CO)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CCCNC(N)=N)C(=O)NC(CCCCN)C(=O)NC(C(C)C)C(=O)NC(CC(C)C)C(=O)NCC(=O)NC(CCC(N)=O)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(C)C(=O)NC(CCCNC(N)=N)C(=O)NC(CCCCN)C(=O)NC(CC(C)C)C(=O)NC(CC(C)C)C(=O)NC(CCC(N)=O)C(=O)NC(CC(O)=O)C(=O)NC(C(C)CC)C(=O)NC(CCSC)C(=O)NC(CO)C(=O)NC(CCCNC(N)=N)C(N)=O", description="Test Smiles-MF 42", outcome="C149H246N44O42S"}
         , TestData {function=testMolForm, input="CC(C)NCCCC1(C(N)=O)C2=CC=CC=C2C2=CC=CC=C12", description="Test Smiles-MF 43", outcome="C20H24N2O"}
         , TestData {function=testMolWt, input="CC[C@H](C)[C@H](NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](CC(=O)O)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](Cc1cnc[nH]1)NC(=O)[C@H](CO)NC(=O)[C@H](CCC(=O)N)NC(=O)[C@@H]1CCCN1C(=O)[C@H](CCCCN)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)CNC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]1CSSC[C@@H]2NC(=O)[C@@H](NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CO)NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CSSC[C@H](NC(=O)[C@H](CCCCN)NC(=O)[C@H](CC(=O)N)NC(=O)CNC(=O)[C@H](CCC(=O)N)NC(=O)CNC2=O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)NCC(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CCC(=O)N)C(=O)N1)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H]1CSSC[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C)[C@@H](C)O)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(C)C)C(=O)N1)C(C)C)C(C)C)[C@@H](C)O)[C@@H](C)O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)N)C(=O)O", description="Test Smiles-MW 44", outcome="69634"}
         {--, TestData {function=testMolWt, input="N=C(\\N)/NCCC[C@H](NC(=O)[C@@H]1CCCN1C(=O)[C@H](N)Cc1ccccc1)C(=O)N1CCC[C@H]1C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)C(=O)O", description="Test Smiles-MW 45", outcome="21802"}
         , TestData {function=testMolWt, input="CCNC(=O)C1CCCN1C(=O)C(CCC/N=C(\\N)/N)NC(=O)C(CC(C)C)NC(=O)C(CC(C)C)NC(=O)C(Cc1ccc(O)cc1)NC(=O)C(CO)NC(=O)C(Cc1c[nH]c2ccccc12)NC(=O)C(Cc1cnc[nH]1)NC(=O)C1CCC(=O)N1", description="Test Smiles-MW 46", outcome="12094"}
         , TestData {function=testMolWt, input="CCC(C)C(NC(=O)C(C)NC(=O)C(CC(=O)O)NC(=O)C(C)NC(=O)C(N)Cc1ccc(O)cc1)C(=O)NC(Cc1ccccc1)C(=O)NC(C(C)O)C(=O)NC(CC(=O)N)C(=O)NC(CO)C(=O)NC(Cc1ccc(O)cc1)C(=O)NC(CCCNC(=N)N)C(=O)NC(CCCCN)C(=O)NC(C(C)C)C(=O)NC(CC(C)C)C(=O)NCC(=O)NC(CCC(=O)N)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(C)C(=O)NC(CCCNC(=N)N)C(=O)NC(CCCCN)C(=O)NC(CC(C)C)C(=O)NC(CC(C)C)C(=O)NC(CCC(=O)N)C(=O)NC(CC(=O)O)C(=O)NC(C(C)CC)C(=O)NC(CCSC)C(=O)NC(CO)C(=O)NC(CCCNC(=N)N)C(=O)N", description="Test Smiles-MW 47", outcome="33578"}
         , TestData {function=testMolWt, input="CC(C)CC(NC(=O)C(COC(C)(C)C)NC(=O)C(Cc1ccc(O)cc1)NC(=O)C(CO)NC(=O)C(Cc1c[nH]c2ccccc12)NC(=O)C(Cc1cnc[nH]1)NC(=O)C1CCC(=O)N1)C(=O)NC(CCC/N=C(\\N)/N)C(=O)N1CCCC1C(=O)NNC(=O)N", description="Test Smiles-MW 48", outcome="12694"}
         , TestData {function=testMolWt, input="CC(C)CC(NC(=O)C(NC(=O)C1CSSCC(N)C(=O)NC(CO)C(=O)NC(CC(=O)N)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(C(C)O)C(=O)N1)C(C)C)C(=O)NCC(=O)NC(CCCCN)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(CCC(=O)N)C(=O)NC(CCC(=O)O)C(=O)NC(CC(C)C)C(=O)NC(Cc1cnc[nH]1)C(=O)NC(CCCCN)C(=O)NC(CC(C)C)C(=O)NC(CCC(=O)N)C(=O)NC(C(C)O)C(=O)NC(Cc1ccc(O)cc1)C(=O)N1CCCC1C(=O)NC(CCCNC(=N)N)C(=O)NC(C(C)O)C(=O)NC(CC(=O)N)C(=O)NC(C(C)O)C(=O)NCC(=O)NC(CO)C(=O)NCC(=O)NC(C(C)O)C(=O)N1CCCC1C(=O)N", description="Test Smiles-MW 49", outcome="34318"}
         , TestData {function=testMolWt, input="CC(C)C[C@@H](NC(=O)[C@H](C)NC(=O)CNC(=O)[C@@H](NC=O)C(C)C)C(=O)N[C@@H](C)C(=O)N[C@H](C(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)NCCO", description="Test Smiles-MW 50", outcome="18822"}
         , TestData {function=testMolWt, input="NC(=O)CCC1NC(=O)C(Cc2ccccc2)NC(=O)C(Cc2ccc(O)cc2)NC(=O)CCSSCC(NC(=O)C(CC(=O)N)NC1=O)C(=O)N1CCCC1C(=O)NC(CCC/N=C(\\N)/N)C(=O)NCC(=O)N", description="Test Smiles-MW 51", outcome="10692"}
         , TestData {function=testMolWt, input="CC(=O)O.CC(C)C[C@H](NC(=O)[C@@H](CCCNC(=O)N)NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CO)NC(=O)[C@@H](Cc1cccnc1)NC(=O)[C@@H](Cc1ccc(Cl)cc1)NC(=O)[C@@H](Cc1cc2ccccc2cc1)NC(=O)C)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@H]1C(=O)N[C@H](C)C(=O)N", description="Test Smiles-MW 52", outcome="14910"}
         , TestData {function=testMolWt, input="NC(=O)C1CSSCCC(=O)NC(CCCC/N=C(\\N)/N)C(=O)NCC(=O)NC(CC(=O)O)C(=O)NC(Cc2c[nH]c3ccccc23)C(=O)N2CCCC2C(=O)N1", description="Test Smiles-MW 53", outcome="8319"}
         , TestData {function=testMolWt, input="CCC(C)C1NC(=O)C(Cc2ccc(O)cc2)NC(=O)C(N)CSSCC(NC(=O)C(CC(=O)N)NC(=O)C(CCC(=O)N)NC1=O)C(=O)N1CCCC1C(=O)NC(CCC/N=C(\\N)/N)C(=O)NCC(=O)N", description="Test Smiles-MW 54", outcome="10502"}
         , TestData {function=testMolWt, input="CCCCCCCCCC(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)NC1C(C)OC(=O)[C@H](CC(=O)c2ccccc2N)NC(=O)[C@H](NC(=O)[C@@H](CO)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)[C@@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCCN)NC(=O)CNC1=O)C(C)CC(=O)O", description="Test Smiles-MW 55", outcome="16206"}
         , TestData {function=testMolWt, input="CCC1NC(=O)C([C@H](O)[C@H](C)C/C=C/C)N(C)C(=O)C(C(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(CC(C)C)N(C)C(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)N(C)C(=O)C(NC(=O)C(CC(C)C)N(C)C(=O)CN(C)C1=O)C(C)C", description="Test Smiles-MW 56", outcome="12026"}
         , TestData {function=testMolWt, input="NCCCCC(NC(=O)C1CCCN1C(=O)C1CSSCC(N)C(=O)NC(Cc2ccccc2)C(=O)NC(Cc2ccccc2)C(=O)NC(CCC(=O)N)C(=O)NC(CC(=O)N)C(=O)N1)C(=O)NCC(=O)N", description="Test Smiles-MW 57", outcome="10402"}
         , TestData {function=testMolWt, input="CCC(C)C1NC(=O)C(Cc2ccc(O)cc2)NC(=O)C(N)CSSCC(NC(=O)C(CC(=O)N)NC(=O)C(NC1=O)C(C)O)C(=O)N1CCCC1C(=O)NC(CC(C)C)C(=O)NCC(=O)N", description="Test Smiles-MW 58", outcome="9801"}
         , TestData {function=testMolWt, input="CC(O)C(CO)NC(=O)C1CSSCC(NC(=O)C(N)Cc2ccccc2)C(=O)NC(Cc2ccccc2)C(=O)NC(Cc2c[nH]c3ccccc23)C(=O)NC(CCCCN)C(=O)NC(C(C)O)C(=O)N1", description="Test Smiles-MW 59", outcome="10192"}
         , TestData {function=testMolWt, input="CC(C)C[C@H](NC(=O)[C@@H](CC(=O)N)NC(=O)[C@H](Cc1ccc(O)cc1)N(C)C(=O)[C@H](CO)NC(=O)[C@@H](Cc1cnccc1)NC(=O)[C@@H](Cc1ccc(Cl)cc1)NC(=O)[C@@H](Cc1cc2c(cccc2)cc1)NC(=O)C)C(=O)N[C@@H](CCCCNC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@H](C)C(=O)N", description="Test Smiles-MW 60", outcome="14160"}
         , TestData {function=testMolWt, input="CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](CC(=O)N)NC(=O)[C@H](CCC(=O)N)NC1=O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N", description="Test Smiles-MW 61", outcome="10071"}
         , TestData {function=testMolWt, input="CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@H](Cc1cnc[nH]1)NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)C)[C@@H](C)O)[C@@H](C)CC)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)N)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](C)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(=O)N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccccc1)C(=O)N", description="Test Smiles-MW 62", outcome="44918"}
         , TestData {function=testMolWt, input="Cc1ncc(COP(=O)(O)O)c(C=O)c1O", description="Test Smiles-MW 63", outcome="2471"}
         , TestData {function=testMolWt, input="[Co+3].[C-]#N.CC(CNC(=O)CC[C@]1(C)[C@@H](CC(=O)N)C2[N-]/C/1=C(/C)\\C\\1=N\\C(=C/C/3=N/C(=C(/C)\\C4=N[C@]2(C)[C@@](C)(CC(=O)N)[C@@H]4CCC(=O)N)/[C@@](C)(CC(=O)N)[C@@H]3CCC(=O)N)\\C(C)(C)[C@@H]1CCC(=O)N)OP(=O)(O)OC1[C@H](CO)O[C@@H](C1O)n1cnc2c1cc(C)c(C)c2", description="Test Smiles-MW 64", outcome="13563"}
         , TestData {function=testMolWt, input="Nc1nc(=O)c2c(NCC(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)N2)[nH]1", description="Test Smiles-MW 65", outcome="4454"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1cnc[nH]1)C(=O)O", description="Test Smiles-MW 66", outcome="1551"}
         , TestData {function=testMolWt, input="C[S+](CC[C@H](N)C(=O)O)C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12", description="Test Smiles-MW 67", outcome="3994"}
         , TestData {function=testMolWt, input="CC(=O)C(=O)O", description="Test Smiles-MW 68", outcome="880"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1ccccc1)C(=O)O", description="Test Smiles-MW 69", outcome="1651"}
         , TestData {function=testMolWt, input="OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12", description="Test Smiles-MW 70", outcome="2443"}
         , TestData {function=testMolWt, input="C[N+](C)(C)CCO", description="Test Smiles-MW 71", outcome="1041"}
         , TestData {function=testMolWt, input="NCCCC[C@H](N)C(=O)O", description="Test Smiles-MW 72", outcome="1461"}
         , TestData {function=testMolWt, input="N[C@@H](CCC/N=C(/N)\\N)C(=O)O", description="Test Smiles-MW 73", outcome="1742"}
         , TestData {function=testMolWt, input="OC[C@H](O)[C@H]1OC(=C(O)C1=O)O", description="Test Smiles-MW 74", outcome="1761"}
         , TestData {function=testMolWt, input="NCCCNCCCCNCCCN", description="Test Smiles-MW 75", outcome="2023"}
         , TestData {function=testMolWt, input="N[C@@H](CC(=O)O)C(=O)O", description="Test Smiles-MW 76", outcome="1331"}
         , TestData {function=testMolWt, input="NCCC[C@H](N)C(=O)O", description="Test Smiles-MW 77", outcome="1321"}
         , TestData {function=testMolWt, input="N[C@@H](CCC(=O)N)C(=O)O", description="Test Smiles-MW 78", outcome="1461"}
         , TestData {function=testMolWt, input="Nc1c2ncn([C@@H]3O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1", description="Test Smiles-MW 79", outcome="3472"}
         , TestData {function=testMolWt, input="CC/C=C\\C/C=C\\C/C=C\\CCCCCCCC(=O)O", description="Test Smiles-MW 80", outcome="2784"}
         , TestData {function=testMolWt, input="N[C@@H](CO)C(=O)O", description="Test Smiles-MW 81", outcome="1050"}
         , TestData {function=testMolWt, input="CSCC[C@H](N)C(=O)O", description="Test Smiles-MW 82", outcome="1492"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1ccc(O)cc1)C(=O)O", description="Test Smiles-MW 83", outcome="1811"}
         , TestData {function=testMolWt, input="C[C@H](CCCC(C)(C)O)[C@H]1CC[C@H]2C(=CC=C3C[C@@H](O)C[C@H](O)C3=C)CCC[C@]12C", description="Test Smiles-MW 84", outcome="4166"}
         , TestData {function=testMolWt, input="C/C(=C\\C=C\\C=C(/C)\\C=C\\C=C(/C)\\C=C\\C1=C(C)CC(O)CC1(C)C)/C=C/C=C(\\C)/C=C/C1C(=CC(O)CC1(C)C)C", description="Test Smiles-MW 85", outcome="5688"}
         , TestData {function=testMolWt, input="N[C@@H](CSSC[C@H](N)C(=O)O)C(=O)O", description="Test Smiles-MW 86", outcome="2403"}
         , TestData {function=testMolWt, input="OC(=O)CCC(=O)O", description="Test Smiles-MW 87", outcome="1180"}
         , TestData {function=testMolWt, input="Cc1cc2c(cc1C)n(C[C@@H](O)[C@@H](O)[C@@H](O)CO)c1nc(=O)[nH]c(=O)c1n2", description="Test Smiles-MW 88", outcome="3763"}
         , TestData {function=testMolWt, input="CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O", description="Test Smiles-MW 89", outcome="2212"}
         , TestData {function=testMolWt, input="N[C@@H](CCC(=O)O)C(=O)O", description="Test Smiles-MW 90", outcome="1471"}
         , TestData {function=testMolWt, input="N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O", description="Test Smiles-MW 91", outcome="3073"}
         , TestData {function=testMolWt, input="CCCC(=O)O[C@H](COC(=O)CC)CO[P@@](=O)(O)OC[C@H](N)C(=O)O", description="Test Smiles-MW 92", outcome="3853"}
         , TestData {function=testMolWt, input="NCC(=O)O", description="Test Smiles-MW 93", outcome="750"}
         , TestData {function=testMolWt, input="C[C@H](CCCC(C)(C)O)[C@H]1CCC2[C@]1(C)CCC/C/2=C\\C=C/1\\C[C@H](O)CCC1=C", description="Test Smiles-MW 94", outcome="4006"}
         , TestData {function=testMolWt, input="Cc1ncc(CO)c(C=O)c1O", description="Test Smiles-MW 95", outcome="1671"}
         , TestData {function=testMolWt, input="N=C(\\N)/N(C)CC(=O)O", description="Test Smiles-MW 96", outcome="1311"}
         , TestData {function=testMolWt, input="CC(C)CC(N)C(=O)O", description="Test Smiles-MW 97", outcome="1311"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", description="Test Smiles-MW 98", outcome="2042"}
         , TestData {function=testMolWt, input="N[C@@H](CS)C(=O)O", description="Test Smiles-MW 99", outcome="1211"}
         , TestData {function=testMolWt, input="Cc1c(CCO)sc[n+]1Cc1cnc(C)nc1N", description="Test Smiles-MW 100", outcome="2653"}
         , TestData {function=testMolWt, input="CC(C)[C@@H](C)/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@]1(C)CCC/C/2=C\\C=C/1\\C[C@@H](O)CCC1=C", description="Test Smiles-MW 101", outcome="3966"}
         , TestData {function=testMolWt, input="CCCCC/C=C\\C/C=C\\C/C=C\\CCCCCCC(=O)O", description="Test Smiles-MW 102", outcome="3064"}
         , TestData {function=testMolWt, input="N[C@@H](CCCNC(=O)N)C(=O)O", description="Test Smiles-MW 103", outcome="1751"}
         , TestData {function=testMolWt, input="C[C@@H](O)[C@H](N)C(=O)O", description="Test Smiles-MW 104", outcome="1191"}
         , TestData {function=testMolWt, input="NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H](CO[P@@](=O)(O)O[P@@](=O)(O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O", description="Test Smiles-MW 105", outcome="6654"}
         , TestData {function=testMolWt, input="Nc1nc(=O)c2nc(CNc3ccc(cc3)C(=O)N[C@@H](CCC(=O)O)C(=O)O)cnc2[nH]1", description="Test Smiles-MW 106", outcome="4413"}
         , TestData {function=testMolWt, input="CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O", description="Test Smiles-MW 107", outcome="3024"}
         , TestData {function=testMolWt, input="C[C@H](N)C(=O)O", description="Test Smiles-MW 108", outcome="890"}
         , TestData {function=testMolWt, input="CC(C)[C@H](N)C(=O)O", description="Test Smiles-MW 109", outcome="1171"}
         , TestData {function=testMolWt, input="CC(=CCO)C=CC=C(C)C=CC1=C(C)CCCC1(C)C", description="Test Smiles-MW 110", outcome="2864"}
         , TestData {function=testMolWt, input="CC(C)CCC[C@@H](C)CCC[C@@H](C)CCC[C@]1(C)CCc2c(O1)c(C)c(C)c(O)c2C", description="Test Smiles-MW 111", outcome="4307"}
         , TestData {function=testMolWt, input="Cc1ncc(CO)c(CO)c1O", description="Test Smiles-MW 112", outcome="1691"}
         , TestData {function=testMolWt, input="OC(=O)CCCC[C@@H]1CCSS1", description="Test Smiles-MW 113", outcome="2063"}
         , TestData {function=testMolWt, input="CC[C@H](C)[C@H](N)C(=O)O", description="Test Smiles-MW 114", outcome="1311"}
         , TestData {function=testMolWt, input="COC(=O)C(Cc1ccccc1)NC(=O)C(N)CC(=O)O", description="Test Smiles-MW 115", outcome="2943"}
         , TestData {function=testMolWt, input="CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2C(=CC=C3C[C@@H](O)CCC3=C)CCC[C@]12C", description="Test Smiles-MW 116", outcome="3846"}
         , TestData {function=testMolWt, input="CC1=CC(=O)c2ccccc2C1=O", description="Test Smiles-MW 117", outcome="1721"}
         , TestData {function=testMolWt, input="Nc1c2ncn([C@@H]3O[C@H](CO[P@@](=O)(O)O[P@@](=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O)c2ncn1", description="Test Smiles-MW 118", outcome="5071"}
         , TestData {function=testMolWt, input="OC(=O)[C@@H]1CCCN1", description="Test Smiles-MW 119", outcome="1151"}
         , TestData {function=testMolWt, input="Nc1c2[nH]cnc2ncn1", description="Test Smiles-MW 120", outcome="1351"}
         , TestData {function=testMolWt, input="N[C@@H](CC(=O)N)C(=O)O", description="Test Smiles-MW 121", outcome="1321"}
         , TestData {function=testMolWt, input="CC[C@H](C)C(=O)O[C@H]1C[C@H](O)C=C2C=C[C@H](C)[C@H](CC[C@@H](O)C[C@@H](O)CC(=O)O)[C@@H]12", description="Test Smiles-MW 122", outcome="4245"}
         , TestData {function=testMolWt, input="COCCCCC(=NOCCN)c1ccc(cc1)C(F)(F)F", description="Test Smiles-MW 123", outcome="3183"}
         , TestData {function=testMolWt, input="CCCCC(=O)N(Cc1ccc(cc1)c1ccccc1c1n[nH]nn1)[C@@H](C(C)C)C(=O)O", description="Test Smiles-MW 124", outcome="4355"}
         , TestData {function=testMolWt, input="CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1[C@H]2CCC[C@H]2C[C@H]1C(=O)O", description="Test Smiles-MW 125", outcome="4165"}
         , TestData {function=testMolWt, input="C[C@@H](Cc1cc(O)c(O)cc1)[C@H](C)Cc1cc(O)c(O)cc1", description="Test Smiles-MW 126", outcome="3023"}
         , TestData {function=testMolWt, input="CC1(C)O[C@@H]2C[C@H]3[C@@H]4C[C@H](F)C5=CC(=O)C=C[C@]5(C)[C@H]4[C@@H](O)C[C@]3(C)[C@@]2(O1)C(=O)CO", description="Test Smiles-MW 127", outcome="4344"}
         , TestData {function=testMolWt, input="NCC(CC(=O)O)c1ccc(Cl)cc1", description="Test Smiles-MW 128", outcome="2136"}
         , TestData {function=testMolWt, input="C[C@H](N)Cc1ccccc1", description="Test Smiles-MW 129", outcome="1352"}
         , TestData {function=testMolWt, input="CSCC[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)CCNC(=O)OCC(C)C)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N", description="Test Smiles-MW 130", outcome="7678"}
         , TestData {function=testMolWt, input="CN1CCCC1c1cnccc1", description="Test Smiles-MW 131", outcome="1622"}
         , TestData {function=testMolWt, input="C[C@@H]1O[C@@]2(CS1)CN1CCC2CC1", description="Test Smiles-MW 132", outcome="1993"}
         , TestData {function=testMolWt, input="OC1N=C(c2ccccc2Cl)c2c(NC1=O)ccc(Cl)c2", description="Test Smiles-MW 133", outcome="3211"}
         , TestData {function=testMolWt, input="COC(=O)CCc1ccc(OCC(O)CNC(C)C)cc1", description="Test Smiles-MW 134", outcome="2953"}
         , TestData {function=testMolWt, input="CC(C)C[C@@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)c1nccnc1)[B](O)O", description="Test Smiles-MW 135", outcome="3842"}
         , TestData {function=testMolWt, input="CCC(O)(C=CCl)C#C", description="Test Smiles-MW 136", outcome="1445"}
         , TestData {function=testMolWt, input="C[C@@](Cc1cc(O)c(O)cc1)(NN)C(=O)O", description="Test Smiles-MW 137", outcome="2262"}
         , TestData {function=testMolWt, input="CC(C)(N)Cc1ccccc1", description="Test Smiles-MW 138", outcome="1492"}
         , TestData {function=testMolWt, input="CC(C)NCCCC1(C(=O)N)c2ccccc2c2ccccc12", description="Test Smiles-MW 139", outcome="3084"}
         , TestData {function=testMolWt, input="COc1cccc(c1)[C@@]1(O)CCCC[C@@H]1CN(C)C", description="Test Smiles-MW 140", outcome="2633"}
         , TestData {function=testMolWt, input="O.Nc1c2ncn([C@@H]3O[C@H](CO)[C@@H](O)[C@@H]3O)c2ncn1", description="Test Smiles-MW 141", outcome="2852"}
         , TestData {function=testMolWt, input="CC(C)NCC(O)COc1ccc(CCOCC2CC2)cc1", description="Test Smiles-MW 142", outcome="3074"}
         , TestData {function=testMolWt, input="OC(Cn1cncn1)(Cn1cncn1)c1c(F)cc(F)cc1", description="Test Smiles-MW 143", outcome="3062"}
         , TestData {function=testMolWt, input="CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(=O)C)[C@@H](N)C1", description="Test Smiles-MW 144", outcome="3124"}
         , TestData {function=testMolWt, input="CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O", description="Test Smiles-MW 145", outcome="7339"}
         , TestData {function=testMolWt, input="O.[Co+2].C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(=O)N)[C@H]2[N-]/C/1=C(/C)\\C\\1=N\\C(=C/C/3=N/C(=C(/C)\\C4=N[C@]2(C)[C@@](C)(CC(=O)N)[C@@H]4CCC(=O)N)/[C@@](C)(CC(=O)N)[C@@H]3CCC(=O)N)\\C(C)(C)[C@@H]1CCC(=O)N)OP(=O)([O-])O[C@@H]1[C@@H](CO)O[C@@H]([C@@H]1O)n1cnc2c1cc(C)c(C)c2", description="Test Smiles-MW 146", outcome="13473"}
         , TestData {function=testMolWt, input="Cn1cnc2c1c(=O)n(C)c(=O)n2C", description="Test Smiles-MW 147", outcome="1941"}
         , TestData {function=testMolWt, input="C[N+](C)(C)CCOC(=O)CCC(=O)OCC[N+](C)(C)C", description="Test Smiles-MW 148", outcome="2903"}
         , TestData {function=testMolWt, input="CCCc1nn(C)c2c1[nH]c(nc2=O)c1c(OCC)ccc(c1)S(=O)(=O)N1CCN(C)CC1", description="Test Smiles-MW 149", outcome="4745"}
         , TestData {function=testMolWt, input="CN(CCOc1ccc(NS(=O)(=O)C)cc1)CCc1ccc(NS(=O)(=O)C)cc1", description="Test Smiles-MW 150", outcome="4415"}
         , TestData {function=testMolWt, input="CCc1c(c(N)nc(N)n1)c1ccc(Cl)cc1", description="Test Smiles-MW 151", outcome="2487"}
         , TestData {function=testMolWt, input="CO[C@H]1[C@@H](C[C@@H]2CN3CCc4c([nH]c5c4ccc(OC)c5)[C@H]3C[C@@H]2[C@@H]1C(=O)OC)OC(=O)c1cc(OC)c(OC)c(OC)c1", description="Test Smiles-MW 152", outcome="6086"}
         , TestData {function=testMolWt, input="CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)C(C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)CN(C)[C@H](C)[C@@H](O)[C@]1(C)O", description="Test Smiles-MW 153", outcome="7489"}
         , TestData {function=testMolWt, input="Clc1ccccc1CN1CCc2c(C1)ccs2", description="Test Smiles-MW 154", outcome="2637"}
         , TestData {function=testMolWt, input="OC(C(=O)OC1CC2CCC(C1)[N+]12CCCC1)(c1ccccc1)c1ccccc1", description="Test Smiles-MW 155", outcome="3925"}
         , TestData {function=testMolWt, input="COc1c(cc(cc1)c1cc2c(cc1)cc(cc2)C(=O)O)C12CC3CC(CC(C3)C1)C2", description="Test Smiles-MW 156", outcome="4125"}
         , TestData {function=testMolWt, input="COc1cc(C(O)CNC(=O)CN)c(OC)cc1", description="Test Smiles-MW 157", outcome="2542"}
         , TestData {function=testMolWt, input="CC(C)(C)S(=O)(=O)C[C@H](Cc1ccccc1)C(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@H](CC1CCCCC1)[C@H](O)[C@H](O)C1CC1", description="Test Smiles-MW 158", outcome="6308"}
         , TestData {function=testMolWt, input="COc1c(OC)c(CS(=O)c2nc3c([nH]2)cc(OC(F)F)cc3)ncc1", description="Test Smiles-MW 159", outcome="3833"}
         , TestData {function=testMolWt, input="CC(C)NC(=O)NS(=O)(=O)c1c(Nc2cccc(C)c2)ccnc1", description="Test Smiles-MW 160", outcome="3484"}
         , TestData {function=testMolWt, input="CN(C)CCCC1(OCc2c1ccc(c2)C#N)c1ccc(F)cc1", description="Test Smiles-MW 161", outcome="3243"}
         , TestData {function=testMolWt, input="CN1CCC[C@@H]1Cc1c[nH]c2c1cc(CCS(=O)(=O)c1ccccc1)cc2", description="Test Smiles-MW 162", outcome="3825"}
         , TestData {function=testMolWt, input="OC(=O)[C@@H]1C/C(=C\\C=[N+]\\2/[C@@H](Cc3cc(O)c(O)cc23)C(=O)[O-])/C=C(N1)C(=O)O", description="Test Smiles-MW 163", outcome="3883"}
         , TestData {function=testMolWt, input="COc1c2n(cc(C(=O)O)c(=O)c2cc(F)c1N1C[C@@H]2CCCN[C@@H]2C1)C1CC1", description="Test Smiles-MW 164", outcome="4014"}
         , TestData {function=testMolWt, input="CC[N+](C)(CC)CCOC(=O)C(O)(C1CCCCC1)c1ccccc1", description="Test Smiles-MW 165", outcome="3485"}
         , TestData {function=testMolWt, input="Cc1c(cccc1O)C(=O)N[C@@H](CSc1ccccc1)[C@H](O)CN1C[C@H]2CCCC[C@H]2C[C@H]1C(=O)NC(C)(C)C", description="Test Smiles-MW 166", outcome="5677"}
         , TestData {function=testMolWt, input="CCC(NC(C)C)C(O)c1cc(O)c(O)cc1", description="Test Smiles-MW 167", outcome="2393"}
         , TestData {function=testMolWt, input="CCC1=C(C)CN(C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC2CCC(C)CC2)C1=O", description="Test Smiles-MW 168", outcome="4906"}
         , TestData {function=testMolWt, input="C[C@H]1C[C@H]2[C@@H]3C[C@H](F)C4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(OC(=O)C)C(=O)COC(=O)C", description="Test Smiles-MW 169", outcome="4945"}
         , TestData {function=testMolWt, input="CC(C)(C)NC(=O)[C@@H]1CN(Cc2cnccc2)CCN1C[C@@H](O)C[C@@H](Cc1ccccc1)C(=O)N[C@@H]1[C@H](O)Cc2ccccc12", description="Test Smiles-MW 170", outcome="6137"}
         , TestData {function=testMolWt, input="O.[Gd+3].CNC(=O)CN(CCN(CCN(CC(=O)[O-])CC(=O)NC)CC(=O)[O-])CC(=O)[O-]", description="Test Smiles-MW 171", outcome="5916"}
         , TestData {function=testMolWt, input="N/C(=N/CC1COC2(CCCCC2)O1)/N", description="Test Smiles-MW 172", outcome="2132"}
         , TestData {function=testMolWt, input="CC[C@H](C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12", description="Test Smiles-MW 173", outcome="4045"}
         , TestData {function=testMolWt, input="FC(F)OC(F)(F)C(F)Cl", description="Test Smiles-MW 174", outcome="1844"}
         , TestData {function=testMolWt, input="CN(C)CCn1nnnc1SCC1=C(N2[C@H](SC1)[C@H](NC(=O)Cc1csc(N)n1)C2=O)C(=O)O", description="Test Smiles-MW 175", outcome="5256"}
         , TestData {function=testMolWt, input="CC(C)C[C@H](CN)CC(=O)O", description="Test Smiles-MW 176", outcome="1592"}
         , TestData {function=testMolWt, input="CN1c2c(cc(Cl)cc2)C(=NC(O)C1=O)c1ccccc1", description="Test Smiles-MW 177", outcome="3007"}
         , TestData {function=testMolWt, input="CN1C(CCl)Nc2cc(Cl)c(cc2S1(=O)=O)S(=O)(=O)N", description="Test Smiles-MW 178", outcome="3602"}
         , TestData {function=testMolWt, input="Nc1cc(O)c(cc1)C(=O)O", description="Test Smiles-MW 179", outcome="1531"}
         , TestData {function=testMolWt, input="CCOc1ccccc1O[C@H]([C@@H]1CNCCO1)c1ccccc1", description="Test Smiles-MW 180", outcome="3133"}
         , TestData {function=testMolWt, input="Cc1c(cc(C#N)c(=O)[nH]1)c1ccncc1", description="Test Smiles-MW 181", outcome="2112"}
         , TestData {function=testMolWt, input="BrCCC(=O)N1CCN(CC1)C(=O)CCBr", description="Test Smiles-MW 182", outcome="3560"}
         , TestData {function=testMolWt, input="CCC(C)C1(CC)C(=O)NC(=O)NC1=O", description="Test Smiles-MW 183", outcome="2122"}
         , TestData {function=testMolWt, input="Cc1c2[nH]c(=O)c3c(nccc3)n(C3CC3)c2ncc1", description="Test Smiles-MW 184", outcome="2662"}
         , TestData {function=testMolWt, input="Clc1cc(Cl)c(CO/N=C(/Cn2ccnc2)\\c2c(Cl)cc(Cl)cc2)cc1", description="Test Smiles-MW 185", outcome="4291"}
         , TestData {function=testMolWt, input="CCC(=O)OCC(=O)[C@@]1(OC(=O)CC)[C@H](C)C[C@H]2[C@@H]3[C@H](Cl)CC4=CC(=O)C=C[C@]4(C)[C@H]3[C@@H](O)C[C@]12C", description="Test Smiles-MW 186", outcome="5210"}
         , TestData {function=testMolWt, input="CC(C)CC1(CC=C)C(=O)NC(=O)NC1=O", description="Test Smiles-MW 187", outcome="2242"}
         , TestData {function=testMolWt, input="Nc1c2ncn([C@H]3C[C@H](O)[C@@H](CO)O3)c2nc(Cl)n1", description="Test Smiles-MW 188", outcome="2856"}
         , TestData {function=testMolWt, input="COc1ccccc1OCC(O)CN1CCN(CC(=O)Nc2c(C)cccc2C)CC1", description="Test Smiles-MW 189", outcome="4275"}
         , TestData {function=testMolWt, input="Nc1cc(C(=O)O)c(O)cc1", description="Test Smiles-MW 190", outcome="1531"}
         , TestData {function=testMolWt, input="CN1[C@@H]2CC[C@@H]1CC(C2)OC(c1ccccc1)c1ccccc1", description="Test Smiles-MW 191", outcome="3074"}
         , TestData {function=testMolWt, input="Clc1c(CCN2CCN(CC2)c2nsc3ccccc23)cc2CC(=O)Nc2c1", description="Test Smiles-MW 192", outcome="4129"}
         , TestData {function=testMolWt, input="CC[C@@H](CO)NC(=O)[C@H]1CN(C)[C@@H]2Cc3cn(C)c4cccc(c34)C2=C1", description="Test Smiles-MW 193", outcome="3534"}
         , TestData {function=testMolWt, input="CCNC(=O)N(CCCN(C)C)C(=O)[C@@H]1C[C@H]2[C@@H](Cc3c[nH]c4cccc2c34)N(CC=C)C1", description="Test Smiles-MW 194", outcome="4516"}
         , TestData {function=testMolWt, input="OC[C@H]1O[C@H](C[C@@H]1O)n1cc(I)c(=O)[nH]c1=O", description="Test Smiles-MW 195", outcome="3540"}
         , TestData {function=testMolWt, input="Nc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1", description="Test Smiles-MW 196", outcome="2483"}
         , TestData {function=testMolWt, input="CC(C)N1CCN(CC1)c1ccc(OC[C@H]2CO[C@@](Cn3cncn3)(O2)c2c(Cl)cc(Cl)cc2)cc1", description="Test Smiles-MW 197", outcome="5324"}
         , TestData {function=testMolWt, input="O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1", description="Test Smiles-MW 198", outcome="2522"}
         , TestData {function=testMolWt, input="C[C@H]1C[C@H]2[C@@H]3CC[C@H](C(=O)C)[C@@]3(C)C[C@H](O)[C@@H]2[C@@]2(C)CCC(=O)C=C12", description="Test Smiles-MW 199", outcome="3444"}
         , TestData {function=testMolWt, input="C[C@@H]1[C@H]2[C@H](O)[C@H]3[C@H](N(C)C)C(=O)/C(=C(\\N)/O)/C(=O)[C@@]3(O)C(=O)C2=C(O)c2c1cccc2O", description="Test Smiles-MW 200", outcome="4444"}
         , TestData {function=testMolWt, input="CCC(=C(CC)c1ccc(O)cc1)c1ccc(O)cc1", description="Test Smiles-MW 201", outcome="2683"}
         , TestData {function=testMolWt, input="CN(C)[C@H]1[C@@H]2C[C@H]3C(=C(O)c4c(O)ccc(Cl)c4[C@@]3(C)O)C(=O)[C@]2(O)C(=O)/C(=C(/N)\\O)/C1=O", description="Test Smiles-MW 202", outcome="4788"}
         , TestData {function=testMolWt, input="Clc1ccccc1C(n1ccnc1)(c1ccccc1)c1ccccc1", description="Test Smiles-MW 203", outcome="3448"}
         , TestData {function=testMolWt, input="[Ca+2].CC(=O)[O-].CC(=O)[O-]", description="Test Smiles-MW 204", outcome="1581"}
         , TestData {function=testMolWt, input="Nc1ccc(cc1)S(=O)(=O)N", description="Test Smiles-MW 205", outcome="1722"}
         , TestData {function=testMolWt, input="NC1CONC1=O", description="Test Smiles-MW 206", outcome="1020"}
         , TestData {function=testMolWt, input="Clc1c(Cl)c2c(NC3=NC(=O)CN3C2)cc1", description="Test Smiles-MW 207", outcome="2560"}
         , TestData {function=testMolWt, input="ClCCNC(=O)N(CCCl)N=O", description="Test Smiles-MW 208", outcome="2140"}
         , TestData {function=testMolWt, input="Cc1noc(NS(=O)(=O)c2ccc(N)cc2)c1C", description="Test Smiles-MW 209", outcome="2673"}
         , TestData {function=testMolWt, input="COCCc1ccc(OCC(O)CNC(C)C)cc1", description="Test Smiles-MW 210", outcome="2673"}
         , TestData {function=testMolWt, input="CCN(C(=O)C=CC)c1ccccc1C", description="Test Smiles-MW 211", outcome="2032"}
         , TestData {function=testMolWt, input="Oc1c(Cc2c(O)oc3ccccc3c2=O)c(=O)c2ccccc2o1", description="Test Smiles-MW 212", outcome="3362"}
         , TestData {function=testMolWt, input="CO/N=C(/C(=O)N[C@H]1[C@H]2SCC(=C(N2C1=O)C(=O)O)CSc1nnnn1C)\\c1csc(N)n1", description="Test Smiles-MW 213", outcome="5115"}
         , TestData {function=testMolWt, input="CCCN(CCC)CCc1c2CC(=O)Nc2ccc1", description="Test Smiles-MW 214", outcome="2603"}
         , TestData {function=testMolWt, input="COc1ccc(cc1)/C(=C(\\c1ccc(OC)cc1)/c1ccc(OC)cc1)/Cl", description="Test Smiles-MW 215", outcome="3808"}
         , TestData {function=testMolWt, input="COC(=O)C1=C(C)NC(=C(C1c1cccc2nonc12)C(=O)OC(C)C)C", description="Test Smiles-MW 216", outcome="3713"}
         , TestData {function=testMolWt, input="CC(=O)Nc1c(I)c(C(=O)O)c(I)c(NC(=O)C)c1I", description="Test Smiles-MW 217", outcome="6139"}
         , TestData {function=testMolWt, input="NCCc1ccn[nH]1", description="Test Smiles-MW 218", outcome="1111"}
         , TestData {function=testMolWt, input="CC1(C)O[C@@H]2CO[C@@]3(COS(=O)(=O)N)OC(C)(C)O[C@H]3[C@@H]2O1", description="Test Smiles-MW 219", outcome="3393"}
         , TestData {function=testMolWt, input="CO[C@]1(NC(=O)CSCC#N)[C@H]2SCC(=C(N2C1=O)C(=O)O)CSc1nnnn1C", description="Test Smiles-MW 220", outcome="4715"}
         , TestData {function=testMolWt, input="CCCc1nc(c(n1Cc1ccc(cc1)c1ccccc1c1n[nH]nn1)C(=O)OCc1c(C)oc(=O)o1)C(C)(C)O", description="Test Smiles-MW 221", outcome="5585"}
         , TestData {function=testMolWt, input="COc1c(Nc2c3ccccc3nc3ccccc23)ccc(NS(=O)(=O)C)c1", description="Test Smiles-MW 222", outcome="3934"}
         , TestData {function=testMolWt, input="Cn1c2c([nH]cn2)c(=O)n(C)c1=O", description="Test Smiles-MW 223", outcome="1801"}
         , TestData {function=testMolWt, input="C[C@@H]1CCN([C@H](C1)C(=O)O)C(=O)[C@H](CCC/N=C(/N)\\N)NS(=O)(=O)c1cccc2c1NC[C@H](C)C2", description="Test Smiles-MW 224", outcome="5086"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)cc2)c(I)c1)C(=O)O", description="Test Smiles-MW 225", outcome="6509"}
         , TestData {function=testMolWt, input="CC(C)N(CCC(C(=O)N)(c1ccccc1)c1ccccn1)C(C)C", description="Test Smiles-MW 226", outcome="3394"}
         , TestData {function=testMolWt, input="CCN(CC)CC(=O)Nc1c(C)cccc1C", description="Test Smiles-MW 227", outcome="2343"}
         , TestData {function=testMolWt, input="NCCC(O)(P(=O)(O)O)P(=O)(O)O", description="Test Smiles-MW 228", outcome="2350"}
         , TestData {function=testMolWt, input="CN1CCC[C@@H]1CCO[C@](C)(c1ccccc1)c1ccc(Cl)cc1", description="Test Smiles-MW 229", outcome="3438"}
         , TestData {function=testMolWt, input="C[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1N[C@H]1C=C(CO)[C@H](O)[C@H](O)[C@H]1O", description="Test Smiles-MW 230", outcome="6456"}
         , TestData {function=testMolWt, input="COc1ccc(cc1)C(CN(C)C)C1(O)CCCCC1", description="Test Smiles-MW 231", outcome="2774"}
         , TestData {function=testMolWt, input="[Na+].C[C@]12CC[C@H]3[C@@H](CCc4c3ccc(OS(=O)(=O)[O-])c4)[C@@H]1CCC2=O", description="Test Smiles-MW 232", outcome="3724"}
         , TestData {function=testMolWt, input="CC(C)OC(=O)CCC/C=C\\C[C@H]1[C@@H](O)C[C@@H](O)[C@@H]1/C=C/[C@@H](O)COc1cccc(c1)C(F)(F)F", description="Test Smiles-MW 233", outcome="5005"}
         , TestData {function=testMolWt, input="CC(=O)OCC(=O)[C@@]12OC3(CCCC3)O[C@@H]1C[C@H]1[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]21C", description="Test Smiles-MW 234", outcome="5025"}
         , TestData {function=testMolWt, input="CNCC[C@@H](Oc1ccccc1C)c1ccccc1", description="Test Smiles-MW 235", outcome="2553"}
         , TestData {function=testMolWt, input="C[C@@H](O)[C@H](NC(=O)[C@@H](C)[C@H](O)[C@@H](C)NC(=O)[C@@H](NC(=O)c1c(C)c(N)nc(n1)[C@H](CC(=O)N)NC[C@H](N)C(=O)N)[C@@H](OC1OC(CO)C(O)C(O)C1OC1OC(CO)C(O)C(OC(=O)N)C1O)c1cnc[nH]1)C(=O)NCCc1nc(cs1)c1nc(cs1)C(=O)NCCC[S+](C)C", description="Test Smiles-MW 236", outcome="14155"}
         , TestData {function=testMolWt, input="OC(=O)CCCc1ccc(cc1)N(CCCl)CCCl", description="Test Smiles-MW 237", outcome="3042"}
         , TestData {function=testMolWt, input="CCOC(=O)c1cncn1C(C)c1ccccc1", description="Test Smiles-MW 238", outcome="2442"}
         , TestData {function=testMolWt, input="CN(Cc1cc2c([nH]c(C)nc2=O)cc1)c1ccc(s1)C(=O)N[C@@H](CCC(=O)O)C(=O)O", description="Test Smiles-MW 239", outcome="4584"}
         , TestData {function=testMolWt, input="CCC12CC(=C)C3C(CCC4=CC(=O)CCC34)C1CC[C@@]2(O)C#C", description="Test Smiles-MW 240", outcome="3244"}
         , TestData {function=testMolWt, input="CN1CC[C@]23[C@H]4Oc5c(O)ccc(C[C@@H]1[C@@H]2C=C[C@@H]4O)c35", description="Test Smiles-MW 241", outcome="2853"}
         , TestData {function=testMolWt, input="CCCN1CCCCC1C(=O)Nc1c(C)cccc1C", description="Test Smiles-MW 242", outcome="2744"}
         , TestData {function=testMolWt, input="CCCCN1CCCCC1C(=O)Nc1c(C)cccc1C", description="Test Smiles-MW 243", outcome="2884"}
         , TestData {function=testMolWt, input="Cc1ccccc1N1CCN(CCc2nnc3CCCCn23)CC1", description="Test Smiles-MW 244", outcome="3254"}
         , TestData {function=testMolWt, input="Nc1nc(=O)c2c([nH]1)n(CCC(CO)CO)cn2", description="Test Smiles-MW 245", outcome="2532"}
         , TestData {function=testMolWt, input="C[C@H](Cn1cnc2c(N)ncnc12)OCP(=O)(O)O", description="Test Smiles-MW 246", outcome="2872"}
         , TestData {function=testMolWt, input="Cc1c(C(=O)N[C@H]2[C@H]3SC(C)(C)[C@@H](N3C2=O)C(=O)O)c(no1)c1c(F)cccc1Cl", description="Test Smiles-MW 247", outcome="4538"}
         , TestData {function=testMolWt, input="NCC1CCC(CC1)C(=O)O", description="Test Smiles-MW 248", outcome="1572"}
         , TestData {function=testMolWt, input="C[C@@H](O)[C@@H]1[C@H]2[C@@H](C)C(=C(N2C1=O)C(=O)O)S[C@@H]1CN[C@@H](C1)C(=O)Nc1cccc(c1)C(=O)O", description="Test Smiles-MW 249", outcome="4755"}
         , TestData {function=testMolWt, input="CC[C@]12CC(=C)[C@H]3[C@@H](CCC4=CCCC[C@H]34)[C@@H]1CC[C@@]2(O)C#C", description="Test Smiles-MW 250", outcome="3104"}
         , TestData {function=testMolWt, input="CO[C@]12[C@H]3N[C@H]3CN1C1=C([C@H]2COC(=O)N)C(=O)C(=C(C)C1=O)N", description="Test Smiles-MW 251", outcome="3343"}
         , TestData {function=testMolWt, input="CCC(C)C1(CC=C)C(=O)N(C)C(=O)N(C)C1=O", description="Test Smiles-MW 252", outcome="2523"}
         , TestData {function=testMolWt, input="Cc1cc2c(cc1C(=C)c1ccc(cc1)C(=O)O)C(C)(C)CCC2(C)C", description="Test Smiles-MW 253", outcome="3484"}
         , TestData {function=testMolWt, input="CCCCCCCN(CC)CCCC(O)c1ccc(NS(=O)(=O)C)cc1", description="Test Smiles-MW 254", outcome="3845"}
         , TestData {function=testMolWt, input="CC[C@]1(O)C[C@H]2CN(C1)CCc1c([nH]c3ccccc13)[C@@](C2)(C(=O)OC)c1c(OC)cc2N(C)[C@@H]3[C@]4(CCN5CC=C[C@](CC)([C@@H]45)[C@@H](O)[C@]3(O)C(=O)N)c2c1", description="Test Smiles-MW 255", outcome="7539"}
         , TestData {function=testMolWt, input="NS(=O)(=O)c1c(Cl)ccc(c1)C1(O)NC(=O)c2ccccc12", description="Test Smiles-MW 256", outcome="3387"}
         , TestData {function=testMolWt, input="CCOc1cc2c(cc1)nc(s2)S(=O)(=O)N", description="Test Smiles-MW 257", outcome="2583"}
         , TestData {function=testMolWt, input="CCCC(C)C1(CC)C(=O)NC(=O)NC1=O", description="Test Smiles-MW 258", outcome="2262"}
         , TestData {function=testMolWt, input="CCCC(CCC)C(=O)O", description="Test Smiles-MW 259", outcome="1442"}
         , TestData {function=testMolWt, input="C[C@@H]1NC(=O)[C@@H](N)CNC(=O)[C@@H](NC(=O)/C(=C\\NC(=O)N)/NC(=O)[C@H](CNC(=O)C[C@@H](N)CCCN)NC1=O)[C@H]1CCN=C(N)N1.NCCC[C@H](N)CC(=O)NC[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H](N)CNC(=O)[C@@H](NC(=O)/C(=C\\NC(=O)N)/NC1=O)[C@H]1CCN=C(N)N1", description="Test Smiles-MW 260", outcome="13214"}
         , TestData {function=testMolWt, input="CN(C)CCc1c[nH]c2c1cc(C[C@@H]1COC(=O)N1)cc2", description="Test Smiles-MW 261", outcome="2873"}
         , TestData {function=testMolWt, input="CC(=O)Nc1ccc(O)cc1", description="Test Smiles-MW 262", outcome="1511"}
         , TestData {function=testMolWt, input="COc1c(OCCCN2CCOCC2)cc2c(Nc3cc(Cl)c(F)cc3)ncnc2c1", description="Test Smiles-MW 263", outcome="4469"}
         , TestData {function=testMolWt, input="COc1c2O[C@H]3[C@@H](O)C=C[C@H]4[C@H]5Cc(cc1)c2[C@@]34CCN5C", description="Test Smiles-MW 264", outcome="2993"}
         , TestData {function=testMolWt, input="CCN1CCN(C(=O)N[C@@H](C(=O)N[C@H]2[C@H]3SC(C)(C)[C@@H](N3C2=O)C(=O)O)c2ccccc2)C(=O)C1=O", description="Test Smiles-MW 265", outcome="5175"}
         , TestData {function=testMolWt, input="CN1C[C@@H](C[C@H]2[C@H]1Cc1c[nH]c3cccc2c13)C(=O)N[C@]1(C)O[C@@]2(O)[C@@H]3CCCN3C(=O)[C@H](Cc3ccccc3)N2C1=O", description="Test Smiles-MW 266", outcome="5836"}
         , TestData {function=testMolWt, input="CN(C)CC/C=C/1\\c2ccccc2CCc2ccccc12", description="Test Smiles-MW 267", outcome="2774"}
         , TestData {function=testMolWt, input="OC[C@H]1O[C@H](C[C@@H]1O)n1cc(F)c(=O)[nH]c1=O", description="Test Smiles-MW 268", outcome="2461"}
         , TestData {function=testMolWt, input="Cc1ccc(cc1)C(=O)c1cc(c(O)c(O)c1)[N+](=O)[O-]", description="Test Smiles-MW 269", outcome="2732"}
         , TestData {function=testMolWt, input="C[C@H]1C[C@H]2[C@@H]3CC[C@](O)(C(=O)C)[C@@]3(C)C[C@H](O)[C@]2(F)[C@@]2(C)C=CC(=O)C=C12", description="Test Smiles-MW 270", outcome="3764"}
         , TestData {function=testMolWt, input="[Fe+4].[N-]=O.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N", description="Test Smiles-MW 271", outcome="2159"}
         , TestData {function=testMolWt, input="[Ca+2].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)C(=O)[O-].OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(O)C(=O)[O-]", description="Test Smiles-MW 272", outcome="4904"}
         , TestData {function=testMolWt, input="CN1CC[C@]23[C@H]4Oc5c(O)ccc(C[C@@H]1[C@@H]2CCC4=O)c35", description="Test Smiles-MW 273", outcome="2853"}
         , TestData {function=testMolWt, input="COc1cc2c(cc1)n(C(=O)c1ccc(Cl)cc1)c(C)c2CC(=O)O", description="Test Smiles-MW 274", outcome="3577"}
         , TestData {function=testMolWt, input="CCC(CO)NCCNC(CC)CO", description="Test Smiles-MW 275", outcome="2043"}
         , TestData {function=testMolWt, input="N=C(\\N=C(\\N)/N)/N(C)C", description="Test Smiles-MW 276", outcome="1291"}
         , TestData {function=testMolWt, input="CC(C)[N+]1(C)[C@@H]2CC[C@@H]1CC(C2)OC(=O)C(CO)c1ccccc1", description="Test Smiles-MW 277", outcome="3324"}
         , TestData {function=testMolWt, input="CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1", description="Test Smiles-MW 278", outcome="3094"}
         , TestData {function=testMolWt, input="CN1CCN(CC1)C1=c2cc(C)sc2=Nc2ccccc2N1", description="Test Smiles-MW 279", outcome="3124"}
         , TestData {function=testMolWt, input="CC(C)NCC(O)COc1ccc(CC(=O)N)cc1", description="Test Smiles-MW 280", outcome="2663"}
         , TestData {function=testMolWt, input="NC(=O)N/N=C/c1ccc(o1)[N+](=O)[O-]", description="Test Smiles-MW 281", outcome="1981"}
         , TestData {function=testMolWt, input="CC[C@@H]1/C=C(/C)\\C[C@H](C)C[C@H](OC)[C@H]2O[C@](O)([C@H](C)C[C@@H]2OC)C(=O)C(=O)N2CCCC[C@H]2C(=O)O[C@@H]([C@@H](C)[C@@H](O)CC1=O)/C(=C/[C@@H]1CC[C@H](Cl)[C@@H](C1)OC)/C", description="Test Smiles-MW 282", outcome="8104"}
         , TestData {function=testMolWt, input="COc1cc2c(cc1)nc([nH]2)S(=O)Cc1ncc(C)c(OC)c1C", description="Test Smiles-MW 283", outcome="3454"}
         , TestData {function=testMolWt, input="NC(=O)c1nccnc1", description="Test Smiles-MW 284", outcome="1231"}
         , TestData {function=testMolWt, input="CN1CCCC(CC2c3ccccc3Sc3ccccc23)C1", description="Test Smiles-MW 285", outcome="3094"}
         , TestData {function=testMolWt, input="OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1", description="Test Smiles-MW 286", outcome="3888"}
         , TestData {function=testMolWt, input="CC(C)(C)c1ccc(cc1)C(O)CCCN1CCC(CC1)C(O)(c1ccccc1)c1ccccc1", description="Test Smiles-MW 287", outcome="4716"}
         , TestData {function=testMolWt, input="COc1ccc(cc1)[C@@H]1Sc2ccccc2N(CCN(C)C)C(=O)[C@@H]1OC(=O)C", description="Test Smiles-MW 288", outcome="4145"}
         , TestData {function=testMolWt, input="CNCCCC1c2ccccc2C=Cc2ccccc12", description="Test Smiles-MW 289", outcome="2633"}
         , TestData {function=testMolWt, input="Nc1ccc(cc1)C(=O)NCC(=O)O", description="Test Smiles-MW 290", outcome="1941"}
         , TestData {function=testMolWt, input="COc1c(OC)cc2c(N)nc(nc2c1)N(C)CCCNC(=O)C1CCCO1", description="Test Smiles-MW 291", outcome="3894"}
         , TestData {function=testMolWt, input="CN1C(=O)OC(C)(C)C1=O", description="Test Smiles-MW 292", outcome="1431"}
         , TestData {function=testMolWt, input="[O-][N+](=O)c1c(ccc(c1)C(F)(F)F)C(=O)C1C(=O)CCCC1=O", description="Test Smiles-MW 293", outcome="3292"}
         , TestData {function=testMolWt, input="CN1c2c(cc(Cl)cc2)N(c2ccccc2)C(=O)CC1=O", description="Test Smiles-MW 294", outcome="3007"}
         , TestData {function=testMolWt, input="N=c\\1/nc(cc(N)n1O)N1CCCCC1", description="Test Smiles-MW 295", outcome="2092"}
         , TestData {function=testMolWt, input="CC(=O)[C@@]1(O)CC[C@H]2[C@@H]3C=C(C)C4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C", description="Test Smiles-MW 296", outcome="3424"}
         , TestData {function=testMolWt, input="Nc1nc(=S)c2c([nH]1)nc[nH]2", description="Test Smiles-MW 297", outcome="1671"}
         , TestData {function=testMolWt, input="CC[C@@H](CO)NC(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(c34)C2=C1", description="Test Smiles-MW 298", outcome="3394"}
         , TestData {function=testMolWt, input="CC(C)(C)c1ccc(CN2CCN(CC2)C(c2ccccc2)c2ccc(Cl)cc2)cc1", description="Test Smiles-MW 299", outcome="4330"}
         , TestData {function=testMolWt, input="C[C@H]1[C@H](NC(=O)/C(=N/OC(C)(C)C(=O)O)/c2csc([NH3+])n2)C(=O)N1S(=O)(=O)[O-]", description="Test Smiles-MW 300", outcome="4354"}
         , TestData {function=testMolWt, input="Clc1cc2c(oc(=O)[nH]2)cc1", description="Test Smiles-MW 301", outcome="1695"}
         , TestData {function=testMolWt, input="CCC1(CCC(=O)NC1=O)c1ccc(N)cc1", description="Test Smiles-MW 302", outcome="2322"}
         , TestData {function=testMolWt, input="OC(C1CCCCN1)c1cc(nc2c1cccc2C(F)(F)F)C(F)(F)F", description="Test Smiles-MW 303", outcome="3783"}
         , TestData {function=testMolWt, input="Nc1ccc(cc1)S(=O)(=O)Nc1ncccn1", description="Test Smiles-MW 304", outcome="2502"}
         , TestData {function=testMolWt, input="CC(O)C(O)C1CNc2c(N1)c(=O)nc(N)[nH]2", description="Test Smiles-MW 305", outcome="2412"}
         , TestData {function=testMolWt, input="CCC1=CC2CN(C1)Cc1c([nH]c3ccccc13)[C@@](C2)(C(=O)OC)c1c(OC)cc2N(C)[C@@H]3[C@]4(CCN5CC=C[C@](CC)([C@@H]45)[C@@H](OC(=O)C)[C@]3(O)C(=O)OC)c2c1", description="Test Smiles-MW 306", outcome="7789"}
         , TestData {function=testMolWt, input="CCCCCOc1ccc(cc1)c1ccc(cc1)c1ccc(cc1)C(=O)N[C@H]1C[C@@H](O)[C@@H](O)NC(=O)[C@@H]2[C@@H](O)[C@@H](C)CN2C(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC1=O)[C@@H](C)O)[C@H](O)[C@@H](O)c1ccc(O)cc1)[C@@H](C)O", description="Test Smiles-MW 307", outcome="11402"}
         , TestData {function=testMolWt, input="CN1CCN(CC1)C1=c2ccccc2=Nc2c(N1)cc(Cl)cc2", description="Test Smiles-MW 308", outcome="3268"}
         , TestData {function=testMolWt, input="O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.O.[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].OS(=O)(=O)OC[C@H]1O[C@@H](O[C@]2(COS(=O)(=O)O)O[C@H](OS(=O)(=O)O)[C@@H](OS(=O)(=O)O)[C@@H]2OS(=O)(=O)O)[C@H](OS(=O)(=O)O)[C@@H](OS(=O)(=O)O)[C@@H]1OS(=O)(=O)O", description="Test Smiles-MW 309", outcome="15719"}
         , TestData {function=testMolWt, input="CC1CN(CCN1)c1c(F)c(C)c2c(=O)c(cn(C3CC3)c2c1)C(=O)O", description="Test Smiles-MW 310", outcome="3593"}
         , TestData {function=testMolWt, input="CN(C)CCOC(C)(c1ccccc1)c1ccccn1", description="Test Smiles-MW 311", outcome="2703"}
         , TestData {function=testMolWt, input="CC[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@H]34)[C@@H]1CC[C@@]2(O)C#C", description="Test Smiles-MW 312", outcome="3124"}
         , TestData {function=testMolWt, input="NC[C@H](O)c1cc(O)c(O)cc1", description="Test Smiles-MW 313", outcome="1691"}
         , TestData {function=testMolWt, input="Nc1nc(=O)n(C[C@@H](CO)OCP(=O)(O)O)cc1", description="Test Smiles-MW 314", outcome="2791"}
         , TestData {function=testMolWt, input="CN1CCN2C(C1)c1ccccc1Cc1c2nccc1", description="Test Smiles-MW 315", outcome="2653"}
         , TestData {function=testMolWt, input="CCCC(C)(COC(=O)N)COC(=O)N", description="Test Smiles-MW 316", outcome="2182"}
         , TestData {function=testMolWt, input="CCSc1cc2c(Sc3ccccc3N2CCCN2CCN(C)CC2)cc1", description="Test Smiles-MW 317", outcome="3996"}
         , TestData {function=testMolWt, input="CC(C)(C)NC[C@H](O)COc1nsnc1N1CCOCC1", description="Test Smiles-MW 318", outcome="3164"}
         , TestData {function=testMolWt, input="CCCCC[C@@H](O)CC[C@@H]1[C@@H](O)C[C@H]2Cc3c(C[C@@H]12)cccc3OCC(=O)O", description="Test Smiles-MW 319", outcome="3905"}
         , TestData {function=testMolWt, input="ClCC1CO1.NCCNCCNCCNCCN", description="Test Smiles-MW 320", outcome="2818"}
         , TestData {function=testMolWt, input="OC(CCN1CCCCC1)(C1CCCCC1)c1ccccc1", description="Test Smiles-MW 321", outcome="3014"}
         , TestData {function=testMolWt, input="O=C1N(C[C@@H]2CCCc3c2c1ccc3)[C@@H]1CN2CCC1CC2", description="Test Smiles-MW 322", outcome="2964"}
         , TestData {function=testMolWt, input="CC(=O)[C@H]1CC[C@H]2[C@@H]3C=CC4=CC(=O)CC[C@@]4(C)[C@@H]3CC[C@]12C", description="Test Smiles-MW 323", outcome="3124"}
         , TestData {function=testMolWt, input="CC(N)COc1c(C)cccc1C", description="Test Smiles-MW 324", outcome="1792"}
         , TestData {function=testMolWt, input="C[C@@H](CN1CC(=O)NC(=O)C1)N1CC(=O)NC(=O)C1", description="Test Smiles-MW 325", outcome="2682"}
         , TestData {function=testMolWt, input="CCOC(=O)C1=C(COCCN)NC(=C(C1c1ccccc1Cl)C(=O)OC)C", description="Test Smiles-MW 326", outcome="4088"}
         , TestData {function=testMolWt, input="Nc1c2CCCCc2nc2ccccc12", description="Test Smiles-MW 327", outcome="1982"}
         , TestData {function=testMolWt, input="CN1CCCN=C1COC(=O)C(O)(C1CCCCC1)c1ccccc1", description="Test Smiles-MW 328", outcome="3444"}
         , TestData {function=testMolWt, input="Nc1nc(N)c2nc(c3ccccc3)c(N)nc2n1", description="Test Smiles-MW 329", outcome="2532"}
         , TestData {function=testMolWt, input="CCCCC(=O)OCC(=O)[C@@]1(O)C[C@H](OC2CC(NC(=O)C(F)(F)F)C(O)C(C)O2)c2c(C1)c(O)c1C(=O)c3c(C(=O)c1c2O)c(OC)ccc3", description="Test Smiles-MW 330", outcome="7236"}
         , TestData {function=testMolWt, input="OC(CCN1CCCC1)(C1CCCCC1)c1ccccc1", description="Test Smiles-MW 331", outcome="2874"}
         , TestData {function=testMolWt, input="CNC[C@H](O)c1cc(O)ccc1", description="Test Smiles-MW 332", outcome="1672"}
         , TestData {function=testMolWt, input="CCOC(=O)n1ccn(C)c1=S", description="Test Smiles-MW 333", outcome="1862"}
         , TestData {function=testMolWt, input="C[C@H]1O[C@H](C[C@H](O)[C@@H]1O)O[C@H]1[C@@H](O)C[C@H](O[C@H]2[C@@H](O)C[C@H](O[C@H]3CC[C@@]4(C)[C@H](CC[C@@H]5[C@@H]4C[C@@H](O)[C@]4(C)C(CC[C@]54O)C4=CC(=O)OC4)C3)O[C@@H]2C)O[C@@H]1C", description="Test Smiles-MW 334", outcome="7809"}
         , TestData {function=testMolWt, input="CCN1CCCC1CNC(=O)c1c(OC)ccc(c1)S(=O)(=O)N", description="Test Smiles-MW 335", outcome="3414"}
         , TestData {function=testMolWt, input="CCN(CC)C(C)CN1c2ccccc2Sc2ccccc12", description="Test Smiles-MW 336", outcome="3124"}
         , TestData {function=testMolWt, input="COCCOC(=O)C1=C(C)NC(=C(C1c1cc(ccc1)[N+](=O)[O-])C(=O)OC(C)C)C", description="Test Smiles-MW 337", outcome="4184"}
         , TestData {function=testMolWt, input="C[C@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(Cl)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO", description="Test Smiles-MW 338", outcome="4089"}
         , TestData {function=testMolWt, input="CCCC(C)(COC(=O)N)COC(=O)NC(C)C", description="Test Smiles-MW 339", outcome="2603"}
         , TestData {function=testMolWt, input="CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C", description="Test Smiles-MW 340", outcome="3144"}
         , TestData {function=testMolWt, input="C[C@@H](N)[C@@H](O)c1ccccc1", description="Test Smiles-MW 341", outcome="1512"}
         , TestData {function=testMolWt, input="CNC(=O)c1nccc(Oc2ccc(NC(=O)Nc3cc(c(Cl)cc3)C(F)(F)F)cc2)c1", description="Test Smiles-MW 342", outcome="4648"}
         , TestData {function=testMolWt, input="OC(Cn1ccnc1)(P(=O)(O)O)P(=O)(O)O", description="Test Smiles-MW 343", outcome="2720"}
         , TestData {function=testMolWt, input="COc1cc(OC)c(Cl)c2c1C(=O)[C@]1(O2)[C@H](C)CC(=O)C=C1OC", description="Test Smiles-MW 344", outcome="3527"}
         , TestData {function=testMolWt, input="COC(=O)C1=C(C)NC(=C(C1c1ccccc1[N+](=O)[O-])C(=O)OCC(C)C)C", description="Test Smiles-MW 345", outcome="3884"}
         , TestData {function=testMolWt, input="CN1CCN(CC1)C(=O)O[C@@H]1N(C(=O)c2nccnc12)c1ncc(Cl)cc1", description="Test Smiles-MW 346", outcome="3888"}
         , TestData {function=testMolWt, input="CCCC[C@@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)CNC(=O)[C@@H](NC(=O)[C@@H](N)Cc1ccc(OS(=O)(=O)O)cc1)[C@@H](C)O)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N", description="Test Smiles-MW 347", outcome="9800"}
         , TestData {function=testMolWt, input="Cc1nnc2CN=C(c3ccccc3)c3c(ccc(Cl)c3)n12", description="Test Smiles-MW 348", outcome="3087"}
         , TestData {function=testMolWt, input="CN(C)CC[C@@H](c1ccc(Br)cc1)c1ccccn1", description="Test Smiles-MW 349", outcome="3192"}
         , TestData {function=testMolWt, input="CN(C)c1ccc(cc1)/C(=C\\1/C=C\\C(=[N+](/C)\\C)\\C=C1)/c1ccc(cc1)N(C)C", description="Test Smiles-MW 350", outcome="3725"}
         , TestData {function=testMolWt, input="CC(=O)NC1C(O)OC(COS(=O)(=O)O)C(OC2OC(C(OC3OC(CO)C(OC4OC(C(O)C(O)C4OS(=O)(=O)O)C(=O)O)C(OS(=O)(=O)O)C3NS(=O)(=O)O)C(O)C2OS(=O)(=O)O)C(=O)O)C1O", description="Test Smiles-MW 351", outcome="11349"}
         , TestData {function=testMolWt, input="CN1CCN(CC1)C1=Nc2ccccc2Oc2c1cc(Cl)cc2", description="Test Smiles-MW 352", outcome="3278"}
         , TestData {function=testMolWt, input="CCN1CCC[C@H]1CNC(=O)c1c(OC)ccc(Br)c1OC", description="Test Smiles-MW 353", outcome="3712"}
         , TestData {function=testMolWt, input="C[C@H](O)[C@H](C)[C@@H]1O[C@H]1C[C@H]1CO[C@@H](C/C(=C/C(=O)OCCCCCCCCC(=O)O)/C)[C@H](O)[C@@H]1O", description="Test Smiles-MW 354", outcome="5006"}
         , TestData {function=testMolWt, input="C[N+](C)(C)CCOC(=O)N", description="Test Smiles-MW 355", outcome="1471"}
         , TestData {function=testMolWt, input="CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1", description="Test Smiles-MW 356", outcome="3574"}
         , TestData {function=testMolWt, input="CCCN[C@@H]1CCc2c(C1)sc(N)n2", description="Test Smiles-MW 357", outcome="2113"}
         , TestData {function=testMolWt, input="CC(=O)c1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1", description="Test Smiles-MW 358", outcome="3243"}
         , TestData {function=testMolWt, input="CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O", description="Test Smiles-MW 359", outcome="3494"}
         , TestData {function=testMolWt, input="[I-].[I-].COc1c2Oc3c(OC)cc4CC[N+](C)(C)[C@@H](Cc5ccc(Oc6c7[C@@H](Cc(cc1)c2)[N+](C)(C)CCc7cc(OC)c6OC)cc5)c4c3", description="Test Smiles-MW 360", outcome="9066"}
         , TestData {function=testMolWt, input="CC1(C)S[C@@H]2[C@H](NC(=O)COc3ccccc3)C(=O)N2[C@H]1C(=O)O", description="Test Smiles-MW 361", outcome="3503"}
         , TestData {function=testMolWt, input="CCC[C@@H](C)C1(CC=C)C(=O)NC(=O)NC1=O", description="Test Smiles-MW 362", outcome="2382"}
         , TestData {function=testMolWt, input="CCCCN1C[C@H](O)[C@@H](O)[C@H](O)[C@H]1CO", description="Test Smiles-MW 363", outcome="2192"}
         , TestData {function=testMolWt, input="CN(C)CCCN1c2ccccc2Sc2ccccc12", description="Test Smiles-MW 364", outcome="2844"}
         , TestData {function=testMolWt, input="CC(=O)S[C@@H]1CC2=CC(=O)CC[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](CC[C@]43CCC(=O)O4)[C@H]12", description="Test Smiles-MW 365", outcome="4165"}
         , TestData {function=testMolWt, input="COC(=O)C(C1CCCCN1)c1ccccc1", description="Test Smiles-MW 366", outcome="2333"}
         , TestData {function=testMolWt, input="COc1ccccc1OCC(O)COC(=O)N", description="Test Smiles-MW 367", outcome="2412"}
         , TestData {function=testMolWt, input="CN1C2CCC1CC(C2)OC(=O)[C@H](CO)c1ccccc1", description="Test Smiles-MW 368", outcome="2893"}
         , TestData {function=testMolWt, input="CN(C)C(=O)Cc1c(nc2ccc(C)cn12)c1ccc(C)cc1", description="Test Smiles-MW 369", outcome="3073"}
         , TestData {function=testMolWt, input="CC(=O)OCC(CCn1cnc2cnc(N)nc12)COC(=O)C", description="Test Smiles-MW 370", outcome="3213"}
         , TestData {function=testMolWt, input="Cc1ccc(cc1)/C(=C\\CN1CCCC1)/c1ccccn1", description="Test Smiles-MW 371", outcome="2783"}
         , TestData {function=testMolWt, input="CN(N=O)C(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O", description="Test Smiles-MW 372", outcome="2652"}
         , TestData {function=testMolWt, input="CCCCC[C@](C)(O)/C=C/[C@H]1[C@H](O)C[C@H](O)[C@@H]1C/C=C/CCCC(=O)O", description="Test Smiles-MW 373", outcome="3685"}
         , TestData {function=testMolWt, input="Cn1nnnc1SCC1=C(N2[C@H](SC1)C(NC(=O)[C@H](NC(=O)c1c[nH]c(C)cc1=O)c1ccc(O)cc1)C2=O)C(=O)O", description="Test Smiles-MW 374", outcome="6126"}
         , TestData {function=testMolWt, input="ClC1C(Cl)C(Cl)C(Cl)C(Cl)C1Cl", description="Test Smiles-MW 375", outcome="2908"}
         , TestData {function=testMolWt, input="OC[C@H]1O[C@H](C[C@@H]1O)n1cc(c(=O)[nH]c1=O)C(F)(F)F", description="Test Smiles-MW 376", outcome="2962"}
         , TestData {function=testMolWt, input="CN1CCN(CCCN2c3ccccc3Sc3c2cc(Cl)cc3)CC1", description="Test Smiles-MW 377", outcome="3739"}
         , TestData {function=testMolWt, input="CN1CC/C(=C/2\\c3ccccc3C=Cc3ccccc23)/CC1", description="Test Smiles-MW 378", outcome="2873"}
         , TestData {function=testMolWt, input="[N]=O", description="Test Smiles-MW 379", outcome="300"}
         , TestData {function=testMolWt, input="NS(=O)(=O)c1cc2c(NC(Cc3ccccc3)NS2(=O)=O)cc1C(F)(F)F", description="Test Smiles-MW 380", outcome="4214"}
         , TestData {function=testMolWt, input="O=c1ncnc2[nH][nH]cc12", description="Test Smiles-MW 381", outcome="1361"}
         , TestData {function=testMolWt, input="CC(C)(O/N=C(\\C(=O)N[C@H]1[C@H]2SCC(=C(N2C1=O)C(=O)[O-])C[n+]1ccccc1)/c1csc(N)n1)C(=O)O", description="Test Smiles-MW 382", outcome="5465"}
         , TestData {function=testMolWt, input="COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", description="Test Smiles-MW 383", outcome="2903"}
         , TestData {function=testMolWt, input="Nc1nc(=O)n(cc1)[C@@H]1O[C@H](CO)[C@@H](O)C1(F)F", description="Test Smiles-MW 384", outcome="2631"}
         , TestData {function=testMolWt, input="Nc1nc(=O)c2c([nH]1)n(cn2)[C@H]1C[C@H](O)[C@@H](CO)C1=C", description="Test Smiles-MW 385", outcome="2772"}
         , TestData {function=testMolWt, input="C[C@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO", description="Test Smiles-MW 386", outcome="3924"}
         , TestData {function=testMolWt, input="COc1cc(cc(OC)c1O)[C@H]1[C@@H]2[C@H](COC2=O)[C@H](O[C@@H]2O[C@@H]3COC(O[C@H]3[C@H](O)[C@H]2O)c2cccs2)c2cc3c(OCO3)cc12", description="Test Smiles-MW 387", outcome="6566"}
         , TestData {function=testMolWt, input="COc1cccc2c1C(=O)c1c(O)c3c(C[C@](O)(C[C@@H]3O[C@H]3C[C@H](N)[C@@H](O)[C@H](C)O3)C(=O)CO)c(O)c1C2=O", description="Test Smiles-MW 388", outcome="5435"}
         , TestData {function=testMolWt, input="OCC(NC(=O)C(Cl)Cl)C(O)c1ccc(cc1)[N+](=O)[O-]", description="Test Smiles-MW 389", outcome="3231"}
         , TestData {function=testMolWt, input="NC(C(=O)NC1C2CCC(=C(N2C1=O)C(=O)O)Cl)c1ccccc1", description="Test Smiles-MW 390", outcome="3497"}
         , TestData {function=testMolWt, input="Cc1c(OCC(F)(F)F)ccnc1CS(=O)c1nc2ccccc2[nH]1", description="Test Smiles-MW 391", outcome="3693"}
         , TestData {function=testMolWt, input="CNCC(O)c1cc(OC(=O)C(C)(C)C)c(OC(=O)C(C)(C)C)cc1", description="Test Smiles-MW 392", outcome="3514"}
         , TestData {function=testMolWt, input="Fc1ccc(cc1)C(=O)CCCN1CCC(=CC1)n1c(=O)[nH]c2ccccc12", description="Test Smiles-MW 393", outcome="3794"}
         , TestData {function=testMolWt, input="N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O", description="Test Smiles-MW 394", outcome="7768"}
         , TestData {function=testMolWt, input="NC[C@@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)[C@H](N)C[C@H](N)[C@H]3O[C@H]3O[C@H](CN)[C@@H](O)[C@H](O)[C@H]3N)[C@@H]2O)[C@H](N)[C@@H](O)[C@@H]1O", description="Test Smiles-MW 395", outcome="6146"}
         , TestData {function=testMolWt, input="CN(C)[C@H]1[C@@H]2C[C@H]3C(=C(O)c4c(O)ccc(Cl)c4[C@@]3(C)O)C(=O)[C@]2(O)C(=O)/C(=C(\\O)/NCO)/C1=O", description="Test Smiles-MW 396", outcome="5089"}
         , TestData {function=testMolWt, input="CCOC(=O)C1(CCN(C)CC1)c1ccccc1", description="Test Smiles-MW 397", outcome="2473"}
         , TestData {function=testMolWt, input="CCOC(=O)N1CC/C(=C\\2/c3c(CCc4c2nccc4)cc(Cl)cc3)/CC1", description="Test Smiles-MW 398", outcome="3828"}
         , TestData {function=testMolWt, input="CC(=O)OCC1=C(N2[C@H](SC1)[C@H](NC(=O)Cc1cccs1)C2=O)C(=O)O", description="Test Smiles-MW 399", outcome="3964"}
         , TestData {function=testMolWt, input="COc1c(OC)cc2c(N)nc(nc2c1)N1CCN(CC1)C(=O)c1ccco1", description="Test Smiles-MW 400", outcome="3834"}
         , TestData {function=testMolWt, input="CN(C)CCCN1c2ccccc2CCc2ccccc12", description="Test Smiles-MW 401", outcome="2804"}
         , TestData {function=testMolWt, input="COc1c(C)c(C)c(/C=C/C(=C/C=C/C(=C\\C(=O)O)/C)/C)c(C)c1", description="Test Smiles-MW 402", outcome="3264"}
         , TestData {function=testMolWt, input="COC(=O)CCc1c2[nH]c(/C=C/3\\N=C(\\C=c/4\\[nH]/c(=C\\C\\5=N\\C(=C/2)\\C(=C5C)CCC(=O)O)/c(C=C)c4C)/C2=CC=C([C@@H](C(=O)OC)[C@@]32C)C(=O)OC)c1C", description="Test Smiles-MW 403", outcome="7187"}
         , TestData {function=testMolWt, input="COc1cc2c(cc1)cc(CCC(=O)C)cc2", description="Test Smiles-MW 404", outcome="2282"}
         , TestData {function=testMolWt, input="C[N+]1(C)C2CC(CC1C1OC21)OC(=O)C(CO)c1ccccc1", description="Test Smiles-MW 405", outcome="3183"}
         , TestData {function=testMolWt, input="CCC1(CC)C(=O)NC(=O)N(C)C1=O", description="Test Smiles-MW 406", outcome="1982"}
         , TestData {function=testMolWt, input="CCCCCCCCCCCCCCOS(=O)(=O)O", description="Test Smiles-MW 407", outcome="2944"}
         , TestData {function=testMolWt, input="OC(=O)C1CCn2c1ccc2C(=O)c1ccccc1", description="Test Smiles-MW 408", outcome="2552"}
         , TestData {function=testMolWt, input="CC(=C)[C@@H]1C2OC(=O)C1[C@]1(O)CC3O[C@@]43C(=O)O[C@H]2[C@]14C.CC(C)(O)[C@@H]1C2OC(=O)C1[C@]1(O)C[C@H]3O[C@@]43C(=O)O[C@H]2[C@]14C", description="Test Smiles-MW 409", outcome="6025"}
       --}  
         , TestData {function=testMolWt, input="C", description="Test Smiles-MW 1", outcome="160"}
         , TestData {function=testMolWt, input="CC", description="Test Smiles-MW 2", outcome="300"}
         , TestData {function=testMolWt, input="CCC", description="Test Smiles-MW 3", outcome="440"}
         , TestData {function=testMolWt, input="C(C)C", description="Test Smiles-MW 4", outcome="440"}
         , TestData {function=testMolWt, input="C(C(C))", description="Test Smiles-MW 5", outcome="440"}
         , TestData {function=testMolWt, input="C(C(C))C", description="Test Smiles-MW 6", outcome="581"}
         , TestData {function=testMolWt, input="C(C(C))(C)", description="Test Smiles-MW 7", outcome="581"}
         , TestData {function=testMolWt, input="C(C(C))CC(C)", description="Test Smiles-MW 8", outcome="861"}
         , TestData {function=testMolWt, input="C(C(C(C)(C)))", description="Test Smiles-MW 9", outcome="721"}
         , TestData {function=testMolWt, input="C(C)(C)(C)C", description="Test Smiles-MW 10", outcome="721"}
         , TestData {function=testMolWt, input="C=C", description="Test Smiles-MW 11", outcome="280"}
         , TestData {function=testMolWt, input="C=C=C", description="Test Smiles-MW 12", outcome="400"}
         , TestData {function=testMolWt, input="C#C", description="Test Smiles-MW 13", outcome="260"}
         , TestData {function=testMolWt, input="C#CC#C", description="Test Smiles-MW 14", outcome="500"}
         , TestData {function=testMolWt, input="C.C", description="Test Smiles-MW 15", outcome="320"}
         , TestData {function=testMolWt, input="C-C", description="Test Smiles-MW 16", outcome="300"}
         , TestData {function=testMolWt, input="C-C-C-C", description="Test Smiles-MW 17", outcome="581"}
         , TestData {function=testMolWt, input="C.C.C.C", description="Test Smiles-MW 18", outcome="641"}
         , TestData {function=testMolWt, input="N", description="Test Smiles-MW 19", outcome="170"}
         , TestData {function=testMolWt, input="CN", description="Test Smiles-MW 20", outcome="310"}
         , TestData {function=testMolWt, input="C#N", description="Test Smiles-MW 21", outcome="270"}
         , TestData {function=testMolWt, input="OCCN", description="Test Smiles-MW 22", outcome="610"}
         , TestData {function=testMolWt, input="CBr", description="Test Smiles-MW 23", outcome="949"}
         , TestData {function=testMolWt, input="CCl", description="Test Smiles-MW 24", outcome="504"}
         , TestData {function=testMolWt, input="CF", description="Test Smiles-MW 25", outcome="340"}
         , TestData {function=testMolWt, input="OS(=O)(=O)O", description="Test Smiles-MW 26", outcome="980"}
         , TestData {function=testMolWt, input="CC(=O)CC(=O)OC", description="Test Smiles-MW 27", outcome="1161"}
         , TestData {function=testMolWt, input="C1CC1", description="Test Smiles-MW 28", outcome="420"}
         , TestData {function=testMolWt, input="C1CCCC1", description="Test Smiles-MW 29", outcome="701"}
         , TestData {function=testMolWt, input="C2CCCC2", description="Test Smiles-MW 30", outcome="701"}
         , TestData {function=testMolWt, input="C3CCC3", description="Test Smiles-MW 31", outcome="561"}
         , TestData {function=testMolWt, input="C1CC=CC1", description="Test Smiles-MW 32", outcome="681"}  
         , TestData {function=testMolWt, input="C1CC1CC", description="Test Smiles-MW 33", outcome="701"}    
         , TestData {function=testMolWt, input="C1.C1", description="Test Smiles-MW 34", outcome="300"}
         , TestData {function=testMolWt, input="C12CCCC1CCCC2", description="Test Smiles-MW 35", outcome="1242"}
         , TestData {function=testMolWt, input="C123CCCC(CC5CC2)(CCC3)(CC5CC1)", description="Test Smiles-MW 36", outcome="2183"} 
         , TestData {function=testMolWt, input="C=1CCCCCC1", description="Test Smiles-MW 37", outcome="961"}  
         , TestData {function=testMolWt, input="c1ccccc1", description="Test Smiles-MW 38", outcome="781"}  
         , TestData {function=testMolWt, input="n1cccc1", description="Test Smiles-MW 39", outcome="660"}
         , TestData {function=testMolWt, input="cc", description="Test Smiles-MW 40", outcome="280"}
         , TestData {function=testMolWt, input="o1cccc1", description="Test Smiles-MW 41", outcome="680"}
         , TestData {function=testMolWt, input="CCC(C)C(NC(=O)C(C)NC(=O)C(CC(O)=O)NC(=O)C(C)NC(=O)C(N)CC1=CC=C(O)C=C1)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(C(C)O)C(=O)NC(CC(N)=O)C(=O)NC(CO)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC(CCCNC(N)=N)C(=O)NC(CCCCN)C(=O)NC(C(C)C)C(=O)NC(CC(C)C)C(=O)NCC(=O)NC(CCC(N)=O)C(=O)NC(CC(C)C)C(=O)NC(CO)C(=O)NC(C)C(=O)NC(CCCNC(N)=N)C(=O)NC(CCCCN)C(=O)NC(CC(C)C)C(=O)NC(CC(C)C)C(=O)NC(CCC(N)=O)C(=O)NC(CC(O)=O)C(=O)NC(C(C)CC)C(=O)NC(CCSC)C(=O)NC(CO)C(=O)NC(CCCNC(N)=N)C(N)=O", description="Test Smiles-MW 42", outcome="33578"}
         , TestData {function=testMolWt, input="CC(C)NCCCC1(C(N)=O)C2=CC=CC=C2C2=CC=CC=C12", description="Test Smiles-MW 43", outcome="3084"}
         
        ]


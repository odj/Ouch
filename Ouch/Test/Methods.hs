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
        ]

performTests :: [TestData] -> String
performTests [] = ""
performTests (x:xs) = output
    where result = (function x) (input x)
          output = case result of
              Left s ->  "\n" ++ description x ++ ":\t" ++ "FAILED\t-with error string:\t" ++ s ++ "\n" ++ performTests xs
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
testMolForm s = molecularFormula $ fillMoleculeValence $ makeMoleculeFromSmiles s

testMolWt::String -> Either String String
testMolWt s = output
    where eitherMolWt =  molecularWeight $ fillMoleculeValence $ makeMoleculeFromSmiles s
          output = case eitherMolWt of
              Left mw   -> Left mw
              Right mw  -> Right (show $ floor (10 * mw))


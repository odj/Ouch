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
     , performTests
     , makeTestFromString
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
import Data.List as List
import Data.Either
import Data.Maybe
import Data.Char as Char


-- Data structure to hold test/result pairs and descriptions.  All 'functions' should
-- evaluate as (String -> Either String String) 
data TestData = TestData { function :: (String -> Either String String)  
                         , description :: String  
                         , input :: String
                         , outcome :: String    
                         } deriving (Show)

makeTestFromString :: String -> TestData
makeTestFromString "" = TestData {function=testTest, description="Empty test", input="", outcome=""}
makeTestFromString s = TestData {function=func, description=l3, input=l2, outcome=l4}
    where (l1:l2:l3:l4:_) = parseAtTab s
          func | l1 == "testTest" = testTest
               | l1 == "testFail" = testFail
               | l1 == "testMolForm" = testMolForm
               | l1 == "testMolWt" = testMolWt
               | otherwise = testTest


performTests :: [TestData] -> (String, String)
performTests [] = ("", "")
performTests td = (summary, errorLog)
   where summary = "\tPassed: " ++ show ((length $ rights results) - (length $ lines errorLog)) ++ "\n"
                   ++ "\tFailed: " ++ show ((length $ lefts results) + (length $ lines errorLog)) ++ "\n"
                   ++ "\n"  ++ "\n++++++++++++++++++++\nPerformed " ++ show (length td) ++ " tests.\n--------------------\n"
         errorLog = detail td results
         results = List.map (\a -> (function a) (input a)) td
         detail [] _ = ""
         detail _ [] = ""
         detail (t:ts) (r:rs) = case r of
             Left s ->  description t ++ ":\t" ++ "FAILED\t-with error string:\t" ++ s ++ "\n" ++ detail ts rs
             Right s -> if s == (outcome t)
                      then detail ts rs
                      else description t ++ ":\t" ++ "FAILED\t-with output: " ++ s ++ " || should get: " ++ (outcome t) ++ "\n" ++ detail ts rs

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

parseAtTab :: String -> [String]
parseAtTab s =  case dropWhile Char.isSpace s of
                     "" -> []
                     s' -> w : parseAtTab s''
                           where (w, s'') = break (=='\t') s'



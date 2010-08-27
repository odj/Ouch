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
     , testArray
     , testAtomCount
     , testTest
     , testFail
     , testMolForm
     ) where

{-# LANGUAGE RecordWildCards, CPP #-}
import Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Property.Composition
import Ouch.Property.Property
import Ouch.Input.Smiles
import Ouch.Data.Atom
import Data.List as List
import Data.Either
import Data.Maybe
import Data.Char as Char
import System.IO
import System.Environment
import Data.Time.Clock


-- Data structure to hold test/result pairs and descriptions.  All 'functions' should
-- evaluate as (String -> Either String String)
data TestData = TestData { function :: (String -> Either String String)
                         , description :: String
                         , input :: String
                         , outcome :: String
                         }


parseAtTab :: String -> [String]
parseAtTab s =  case dropWhile Char.isSpace s of
                     "" -> []
                     s' -> w : parseAtTab s''
                           where (w, s'') = break (=='\t') s'



makeTestFromString :: String -> TestData
makeTestFromString "" = TestData {function=testTest, description="Empty test", input="", outcome=""}
makeTestFromString s = TestData {function=func, description=l3, input=l2, outcome=l4}
    where (l1:l2:l3:l4:_) = parseAtTab s
          func | l1 == "testTest"    = testTest
               | l1 == "testFail"    = testFail
               | l1 == "testMolForm" = testMolForm
               | l1 == "testMolWt"   = testMolWt
               | l1 == "testAtomCount" = testAtomCount
               | l1 == "testHeavyCount" = testHeavyCount
               | otherwise = testTest


performTests :: [TestData] -> (String, String)
performTests [] = ("", "")
performTests td = (summary, errorLog)
   where summary = "++++++++++++++++++++++\nPerforming "
                   ++ show (length td)
                   ++ " tests.\n----------------------\n"
                   ++ "\tPassed: " ++ show ((length td) - (length $ lines errorLog)) ++ "\n"
                   ++ "\tFailed: " ++ show ((length $ lines errorLog))
         errorLog = detail td results
         results = List.map (\a -> (function a) (input a)) td
         detail [] _ = ""
         detail _ [] = ""
         detail (t:ts) (r:rs) = case r of
             Left s ->  description t ++ ":\t" ++ "FAILED\t-with error string:\t"
                                               ++ s ++ "\n" ++ detail ts rs
             Right s -> if s == (outcome t)
                      then detail ts rs
                      else description t ++ ":\t" ++ "FAILED\t-with output: "
                                                  ++ s ++ " || should get: " ++ (outcome t)
                                                  ++ "\n" ++ detail ts rs

testArray :: [(String, [TestData])] -> IO ()
testArray []   = return ()
testArray (x:xs) = do
  let (summary, errorLog) = performTests $ snd x
  time1 <- getCurrentTime
  putStrLn $ "\n\nTest file: " ++ fst x
  putStrLn $ performTests (snd x) `seq` summary
  time2 <- getCurrentTime
  putStrLn $ "\t" ++ (show $ diffUTCTime time2 time1)
                  ++ " seconds."
  appendFile "errorLog.txt" $  "\nTest file: " ++ (fst x) ++ "\n"
                               ++ "=====================================\n"
                               ++ errorLog
  testArray xs



-- Simple test function
testTest::String -> Either String String
testTest s = Right s

-- Simple test fail function
testFail::String -> Either String String
testFail s = Left s

-- Test smiles to formula
testMolForm :: String -> Either String String
testMolForm s = case molecularFormula $ makeMoleculeFromSmiles s of
  Nothing   -> Left "Unable to generate Molecular Formula Property"
  Just prop -> Right mf
    where mf = case (value prop) of StringValue str -> str

testMolWt :: String -> Either String String
testMolWt s = case molecularWeight $ makeMoleculeFromSmiles s of
  Nothing   -> Left "Unable to generate Molecular Formula Property"
  Just prop -> Right $ show $ floor (10 * mw)
    where mw = case (value prop) of DoubleValue num -> num

testAtomCount :: String -> Either String String
testAtomCount s = case atomCount $ makeMoleculeFromSmiles s of
  Nothing   -> Left "Unable to generate Molecular Formula Property"
  Just prop -> Right ac
    where ac = case (value prop) of IntegerValue i -> show i

testHeavyCount :: String -> Either String String
testHeavyCount s = case heavyAtomCount $ makeMoleculeFromSmiles s of
  Nothing   -> Left "Unable to generate Molecular Formula Property"
  Just prop -> Right ac
    where ac = case (value prop) of IntegerValue i -> show i






{-------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Module      :  Ouch.Test.Methods
--  Maintainer  :  Orion Jankowski
--  Stability   :  Unstable
--  Portability :


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

--------------------------------------------------------------------------------
-------------------------------------------------------------------------------}

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
     , testForm
     , testRoundTrip
     ) where

import Ouch (version)
import Ouch.Structure.Molecule
import Ouch.Structure.Atom
import Ouch.Structure.Bond
import Ouch.Structure.Marker
import Ouch.Property.Composition
import Ouch.Property.Builder
import Ouch.Enumerate.Method
import Ouch.Input.Smiles
import Ouch.Data.Atom
import Ouch.Enumerate.Formula
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
               | l1 == "testEnum" = testEnum
               | l1 == "testForm" = testForm
               | l1 == "testRoundTrip" = testRoundTrip
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
  putStrLn  "\n"
  putStrLn  version
  putStrLn $ "Test file: " ++ fst x
  putStrLn $ performTests (snd x) `seq` summary
  time2 <- getCurrentTime
  putStrLn $ "\t" ++ (show $ diffUTCTime time2 time1)
                  ++ " seconds."
  appendFile "errorLog.txt" $  "\n\n\nTest file: " ++ (fst x) ++ "\n"
                               ++ "=====================================\n"
                               ++ errorLog
  testArray xs



-- Simple test function
testTest :: String -> Either String String
testTest s = Right s

-- Simple test fail function
testFail :: String -> Either String String
testFail s = Left s

-- Little utility function, unsafe for general use.
right p = case p of Right r -> r

-- Test smiles to formula
testMolForm :: String -> Either String String
testMolForm s = Right $ show $ (right $ value molecularFormula) $ readSmi s

testMolWt :: String -> Either String String
testMolWt s = Right $ show $ (right $ value molecularWeight) $ readSmi s

testAtomCount :: String -> Either String String
testAtomCount s = Right $ show $ (right $ value atomCount) $ readSmi s

testHeavyCount :: String -> Either String String
testHeavyCount s = Right $ show $ (right $ value heavyAtomCount) $ readSmi s

testEnum :: String -> Either String String
testEnum s = Right $ show $ length $ [mol] >#> method
                                           >#> method
                                           >#> method
                                           >#> method
                                           >#> method
  where mol = makeScaffoldFromSmiles s
        scaffoldList = map makeScaffoldFromSmiles ["Cl", "OC(=O)C", "C(=O)C", "C1CCCCC1"]
        bondList = replicate 6 Single
        list = zip bondList scaffoldList
        method = Just $ AddMethod Nothing Nothing (\_ _ -> True) list

testForm :: String -> Either String String
testForm s = Right $ show $ List.length $ expand (read s :: Formula)

testRoundTrip :: String -> Either String String
testRoundTrip s = let
  mol1 = read s :: Molecule
  smi1 = show mol1
  mol2 = read smi1 :: Molecule
  smi2 = show mol2
  in Right $  show $ (smi1 == smi2)

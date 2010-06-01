module Ouch.Data.Atom (
      atomicWeights
    , atomicSymbols
    , atomicNames
    , molecularFormulaElements
    ) where

import qualified Data.Map as M


-- atomicWeights
-- Reference: Pure Appl. Chem.,  81,  2131-2156 (2009) 
atomicWeights = M.fromList [
                 (1, 1.00794)
               , (2, 4.002602)
               , (3, 6.941)
               , (4, 9.012182)
               , (5, 10.811)
               , (6, 12.0107)
               , (7, 14.0067)
               , (8, 15.9994)
               , (9, 18.9984032)
               , (10, 20.1797)
               , (11, 22.98976928)
               , (12, 24.3050)
               , (13, 26.9815386)
               , (14, 28.0855)
               , (15, 30.973762)
               , (16, 32.065)
               , (17, 35.453)
               , (18, 39.948)
               , (19, 39.0983)
               , (20, 40.078)
               , (21, 44.955912)
               , (22, 47.867)
               , (23, 50.9415)
               , (24, 51.9961)
               , (25, 54.938045)
               , (26, 55.845)
               , (27, 58.933195)
               , (28, 58.6934)
               , (29, 63.546)
               , (30, 65.38)
               , (31, 69.723)
               , (32, 72.64)
               , (33, 74.92160)
               , (34, 78.96)
               , (35, 79.904)
               , (36, 83.798)
               , (37, 85.4678)
               , (38, 87.62)
               , (39, 88.90585)
               , (40, 91.224)
               , (41, 92.90638)
               , (42, 95.96)
               , (43, 98)
               , (44, 101.07)
               , (45, 102.90550)
               , (46, 106.42)
               , (47, 107.8682)
               , (48, 112.411)
               , (49, 114.818)
               , (50, 118.710)
               , (51, 121.760)
               , (52, 127.60)
               , (53, 126.90447)
               , (54, 131.293)
               , (55, 132.9054519)
               , (56, 137.327)
               , (57, 138.90547)
               , (58, 140.116)
               , (59, 140.90765)
               , (60, 144.242)
               , (61, 145)
               , (62, 150.36)
               , (63, 151.964)
               , (64, 157.25)
               , (65, 158.92535)
               , (66, 162.500)
               , (67, 164.93032)
               , (68, 167.259)
               , (69, 168.93421)
               , (70, 173.054)
               , (71, 174.9668)
               , (72, 178.49)
               , (73, 180.94788)
               , (74, 183.84)
               , (75, 186.207)
               , (76, 190.23)
               , (77, 192.217)
               , (78, 195.084)
               , (79, 196.966569)
               , (80, 200.59)
               , (81, 204.3833)
               , (82, 207.2)
               , (83, 208.98040)
               , (84, 209)
               , (85, 210)
               , (86, 222)
               , (87, 223)
               , (88, 226)
               , (89, 227)
               , (90, 232.03806)
               , (91, 231.03588)
               , (92, 238.02891)
               , (93, 237)
               , (94, 244)
               , (95, 243)
               , (96, 247)
               , (97, 247)
               , (98, 252)
               , (99, 252)
               , (100, 257)
               , (101, 258)
               , (102, 259)
               , (103, 262)
               , (104, 265)
               , (105, 268)
               , (106, 271)
               , (107, 272)
               , (108, 270)
               , (109, 276)
               , (110, 281)
               , (111, 280)
               , (112, 285)
               , (113, 284)
               , (114, 289)
               , (115, 288)
               , (116, 293)
               , (118, 294)
             ]


-- atomicWeights
-- Reference: Pure Appl. Chem.,  81,  2131-2156 (2009)
atomicSymbols = M.fromList [
              (1, "H")
            , (2, "He")
            , (3, "Li")
            , (4, "Be")
            , (5, "B")
            , (6, "C")
            , (7, "N")
            , (8, "O")
            , (9, "F")
            , (10, "Ne")
            , (11, "Na")
            , (12, "Mg")
            , (13, "Al")
            , (14, "Si")
            , (15, "P")
            , (16, "S")
            , (17, "Cl")
            , (18, "Ar")
            , (19, "K")
            , (20, "Ca")
            , (21, "Sc")
            , (22, "Ti")
            , (23, "V")
            , (24, "Cr")
            , (25, "Mn")
            , (26, "Fe")
            , (27, "Co")
            , (28, "Ni")
            , (29, "Cu")
            , (30, "Zn")
            , (31, "Ga")
            , (32, "Ge")
            , (33, "As")
            , (34, "Se")
            , (35, "Br")
            , (36, "Kr")
            , (37, "Rb")
            , (38, "Sr")
            , (39, "Y")
            , (40, "Zr")
            , (41, "Nb")
            , (42, "Mo")
            , (43, "Tc")
            , (44, "Ru")
            , (45, "Rh")
            , (46, "Pd")
            , (47, "Ag")
            , (48, "Cd")
            , (49, "In")
            , (50, "Sn")
            , (51, "Sb")
            , (52, "Te")
            , (53, "I")
            , (54, "Xe")
            , (55, "Cs")
            , (56, "Ba")
            , (57, "La")
            , (58, "Ce")
            , (59, "Pr")
            , (60, "Nd")
            , (61, "Pm")
            , (62, "Sm")
            , (63, "Eu")
            , (64, "Gd")
            , (65, "Tb")
            , (66, "Dy")
            , (67, "Ho")
            , (68, "Er")
            , (69, "Tm")
            , (70, "Yb")
            , (71, "Lu")
            , (72, "Hf")
            , (73, "Ta")
            , (74, "W")
            , (75, "Re")
            , (76, "Os")
            , (77, "Ir")
            , (78, "Pt")
            , (79, "Au")
            , (80, "Hg")
            , (81, "Tl")
            , (82, "Pb")
            , (83, "Bi")
            , (84, "Po")
            , (85, "At")
            , (86, "Rn")
            , (87, "Fr")
            , (88, "Ra")
            , (89, "Ac")
            , (90, "Th")
            , (91, "Pa")
            , (92, "U")
            , (93, "Np")
            , (94, "Pu")
            , (95, "Am")
            , (96, "Cm")
            , (97, "Bk")
            , (98, "Cf")
            , (99, "Es")
            , (100, "Fm")
            , (101, "Md")
            , (102, "No")
            , (103, "Lr")
            , (104, "Rf")
            , (105, "Db")
            , (106, "Sg")
            , (107, "Bh")
            , (108, "Hs")
            , (109, "Mt")
            , (110, "Ds")
            , (111, "Rg")
            , (112, "Cn")
            , (113, "Uut")
            , (114, "Uuq")
            , (115, "Uup")
            , (116, "Uuh")
            , (118, "Uuo")
        ]
        
-- atomicWeights
-- Reference: Pure Appl. Chem.,  81,  2131-2156 (2009)
atomicNames = M.fromList [
              (1, "Hydrogen")
            , (2, "Helium")
            , (3, "Lithium")
            , (4, "Beryllium")
            , (5, "Boron")
            , (6, "Carbon")
            , (7, "Nitrogen")
            , (8, "Oxygen")
            , (9, "Fluorine")
            , (10, "Neon")
            , (11, "Sodium")
            , (12, "Magnesium")
            , (13, "Aluminium")
            , (14, "Silicon")
            , (15, "Phosphorus")
            , (16, "Sulfur")
            , (17, "Chlorine")
            , (18, "Argon")
            , (19, "Potassium")
            , (20, "Calcium")
            , (21, "Scandium")
            , (22, "Titanium")
            , (23, "Vanadium")
            , (24, "Chromium")
            , (25, "Manganese")
            , (26, "Iron")
            , (27, "Cobalt")
            , (28, "Nickel")
            , (29, "Copper")
            , (30, "Zinc")
            , (31, "Gallium")
            , (32, "Germanium")
            , (33, "Arsenic")
            , (34, "Selenium")
            , (35, "Bromine")
            , (36, "Krypton")
            , (37, "Rubidium")
            , (38, "Strontium")
            , (39, "Yttrium")
            , (40, "Zirconium")
            , (41, "Niobium")
            , (42, "Molybdenum")
            , (43, "Technetium")
            , (44, "Ruthenium")
            , (45, "Rhodium")
            , (46, "Palladium")
            , (47, "Silver")
            , (48, "Cadmium")
            , (49, "Indium")
            , (50, "Tin")
            , (51, "Antimony")
            , (52, "Tellurium")
            , (53, "Iodine")
            , (54, "Xenon")
            , (55, "Caesium")
            , (56, "Barium")
            , (57, "Lanthanum")
            , (58, "Cerium")
            , (59, "Praseodymium")
            , (60, "Neodymium")
            , (61, "Promethium")
            , (62, "Samarium")
            , (63, "Europium")
            , (64, "Gadolinium")
            , (65, "Terbium")
            , (66, "Dysprosium")
            , (67, "Holmium")
            , (68, "Erbium")
            , (69, "Thulium")
            , (70, "Ytterbium")
            , (71, "Lutetium")
            , (72, "Hafnium")
            , (73, "Tantalum")
            , (74, "Tungsten")
            , (75, "Rhenium")
            , (76, "Osmium")
            , (77, "Iridium")
            , (78, "Platinum")
            , (79, "Gold")
            , (80, "Mercury")
            , (81, "Thallium")
            , (82, "Lead")
            , (83, "Bismuth")
            , (84, "Polonium")
            , (85, "Astatine")
            , (86, "Radon")
            , (87, "Francium")
            , (88, "Radium")
            , (89, "Actinium")
            , (90, "Thorium")
            , (91, "Protactinium")
            , (92, "Uranium")
            , (93, "Neptunium")
            , (94, "Plutonium")
            , (95, "Americium")
            , (96, "Curium")
            , (97, "Berkelium")
            , (98, "Californium")
            , (99, "Einsteinium")
            , (100, "Fermium")
            , (101, "Mendelevium")
            , (102, "Nobelium")
            , (103, "Lawrencium")
            , (104, "Rutherfordium")
            , (105, "Dubnium")
            , (106, "Seaborgium")
            , (107, "Bohrium")
            , (108, "Hassium")
            , (109, "Meitnerium")
            , (110, "Darmstadtium")
            , (111, "Roentgenium")
            , (112, "Copernicium")
            , (113, "Ununtrium")
            , (114, "Ununquadium")
            , (115, "Ununpentium")
            , (116, "Ununhexium")
            , (118, "Ununoctium")
            ]
 
molecularFormulaElements = [
        "C"
      , "H"    
      , "N"
      , "O"
      , "Ac"
      , "Al"
      , "Am"
      , "Sb"
      , "Ar"
      , "As"
      , "At"
      , "Ba"
      , "Bk"
      , "Be"
      , "Bi"
      , "Bh"
      , "B"
      , "Br"
      , "Cd"
      , "Cs"
      , "Ca"
      , "Cf"
      , "Ce"
      , "Cl"
      , "Cr"
      , "Co"
      , "Cn"
      , "Cu"
      , "Cm"
      , "Ds"
      , "Db"
      , "Dy"
      , "Es"
      , "Er"
      , "Eu"
      , "Fm"
      , "F"
      , "Fr"
      , "Gd"
      , "Ga"
      , "Ge"
      , "Au"
      , "Hf"
      , "Hs"
      , "He"
      , "Ho"
      , "In"
      , "I"
      , "Ir"
      , "Fe"
      , "Kr"
      , "La"
      , "Lr"
      , "Pb"
      , "Li"
      , "Lu"
      , "Mg"
      , "Mn"
      , "Mt"
      , "Md"
      , "Hg"
      , "Mo"
      , "Nd"
      , "Ne"
      , "Np"
      , "Ni"
      , "Nb"
      , "No"
      , "Os"
      , "Pd"
      , "P"
      , "Pt"
      , "Pu"
      , "Po"
      , "K"
      , "Pr"
      , "Pm"
      , "Pa"
      , "Ra"
      , "Rn"
      , "Re"
      , "Rh"
      , "Rg"
      , "Rb"
      , "Ru"
      , "Rf"
      , "Sm"
      , "Sc"
      , "Sg"
      , "Se"
      , "Si"
      , "Ag"
      , "Na"
      , "Sr"
      , "S"
      , "Ta"
      , "Tc"
      , "Te"
      , "Tb"
      , "Tl"
      , "Th"
      , "Tm"
      , "Sn"
      , "Ti"
      , "W"
      , "Uuh"
      , "Uuo"
      , "Uup"
      , "Uuq"
      , "Uut"
      , "U"
      , "V"
      , "Xe"
      , "Yb"
      , "Y"
      , "Zn"
      , "Zr"
      , "Zr"
      ]
                     
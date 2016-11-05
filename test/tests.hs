-- file: test.hs
--
module Main ( main) where

import System.Environment
import System.Exit
import Text.Printf

import Input ( hSlurp, dataFiles, hSlurpAll )
import Types ( M, V, XMeas (..), HMeas (..), PMeas (..), Prong (..), VHMeas (..)
             , showXMeas, showPMeas, showQMeas, showMMeas, showHMeas )
import Coeff ( w2pt, h2p, h2q, q2p, invMass )
import Matrix ( inv, sw, fromList, fromList2 )
import Fit ( fit )

main :: IO ()
main = getArgs >>= parse

parse ["-h"] = usage   >> exit
parse ["-v"] = version >> exit
parse []     = test ["1"]
parse args   = test args

usage   = putStrLn "Usage: fvt [-vh] [test# ..]"
version = putStrLn "Haskell fvt 0.1"
exit    = exitWith ExitSuccess

-- filter list of tracks of helices etc given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
hFilter :: ( Eq a, Enum a, Num a ) => [b] -> [a] -> [b]
hFilter hl rng =
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

thisFile  = "dat/tr05129e001412.dat"
otherFile = "dat/tr05158e004656.dat"
thirdFile = "dat/tr05166e001984.dat"

showP :: (HMeas -> IO ())
showP = (showPMeas "px,py,pz,E ->" . h2p)
showQ :: (HMeas -> IO ())
showQ =  (showQMeas "pt,pz,𝜑 ,E ->" . h2q)
showMomentum :: HMeas -> IO ()
showMomentum h = do
--  showP h
  showQ h

test :: [String] -> IO ()
test arg =
  case arg of
    ["1"] -> do
          VHMeas v hl <- hSlurp thisFile
          mapM_ showMomentum hl
          let l5 = [0,2,3,4,5] -- these are the tracks supposedly from the tau
          doFitTest v hl l5

    ["2"] -> do
          VHMeas v hl <- hSlurp otherFile
          mapM_ showMomentum hl
          let l5 = [0,1,2,4,5]
          doFitTest v hl l5

-- slurp in all event data files from ./dat and append helices
    ["3"] -> do
          ps <- dataFiles "dat"
          VHMeas v hl <- hSlurpAll ps
          mapM_ showMomentum hl
          doFitTest v hl [0..]
 --         putStrLn $ showXMeas "ok, let's check it" v

doFitTest :: XMeas -> [HMeas] -> [Int] -> IO ()
doFitTest v hl l5 = do
  let showLen xs = show $ length xs
  let showQChi2 (qm, chi2, i) = showQMeas ((printf "q%d chi2 ->%8.1f " (i::Int) (chi2::Double)) ++ "pt,pz,fi,E ->") qm

  putStrLn $ showXMeas "initial vertex position ->" v

  let pl              = map h2p hl
  showMMeas ("Inv Mass " ++ showLen pl ++ " helix") $ invMass pl
  let pl5             = map h2p $ hFilter hl l5
  showMMeas ("Inv Mass " ++ showLen pl5 ++ " helix") $ invMass pl5

  putStrLn            "Fitting Vertex --------------------"
  let Prong n vf ql cl = fit v hl
  putStrLn $ showXMeas "Fitted vertex ->" vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  showMMeas ("Inv Mass " ++ showLen ql ++ " fit") $ invMass $map q2p ql
  let pl5              = map q2p $ hFilter ql l5
  showMMeas ("Inv Mass " ++ showLen pl5 ++ " fit") $ invMass pl5

  putStrLn            "Refitting Vertex-----------------"
  let Prong n vf ql cl = fit v $ hFilter hl l5
  putStrLn $ showXMeas "Refitted vertex ->" vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  showMMeas ("Inv Mass " ++ showLen ql ++ " refit")  $ invMass $ map q2p ql


-- file: test.hs
--
module Main ( main) where

import System.Environment
import System.Exit
import Text.Printf
-- :set -XQuasiQuotes
-- import Data.String.Interpolate

import Input ( hSlurp, dataFiles, hSlurpAll )
import Types (  XMeas (..), HMeas (..), Prong (..), VHMeas (..)
              , showXMeas, showPMeas, showQMeas, showMMeas, showHMeas
              , showXMDist, origin
              , h2p, h2q, q2p
             )
import Coeff ( invMass )
import Fit ( fit, fit' )

main :: IO ()
main = getArgs >>= parse

parse :: [String] -> IO ()
parse ["-h"] = usage   >> exit
parse ["-v"] = version >> exit
parse []     = test ["1"]
parse args   = test args

usage :: IO ()
usage   = putStrLn "Usage: fvt [-vh] [test# ..]"
version :: IO ()
version = putStrLn "Haskell fvt 0.1"
exit :: IO ()
exit    = exitSuccess

-- filter list of tracks of helices etc given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
hFilter :: ( Eq a, Enum a, Num a ) => [b] -> [a] -> [b]
hFilter hl rng =
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

thisFile  = "dat/tr05129e001412.dat"
otherFile = "dat/tr05158e004656.dat"
thirdFile = "dat/tr05166e001984.dat"
mc0File   = "dat/tr00101e007076.dat"
mc1File   = "dat/tr00101e008340.dat"
mc2File   = "dat/tr00101e012040.dat"
cmsFile   = "dat/tandv.dat"

showP :: (HMeas -> IO ())
showP = showPMeas "px,py,pz,E ->" . h2p
showQ :: (HMeas -> IO ())
showQ = showQMeas "pt,pz,fi,E ->" . h2q
showH :: (HMeas -> String)
showH = showHMeas "Helix ->"
showMomentum :: HMeas -> IO ()
showMomentum = showQ
showHelix :: HMeas -> IO ()
showHelix = putStrLn . showH

showProng :: Prong -> IO ()
showProng (Prong _ v ql cl) = do
  let
      showCl :: String -> [Double] -> String
      showCl = foldl (\s x -> s++printf "%8.1g" (x::Double))
      st = showXMDist (printf "chi2tot ->%8.1f, r ->" (sum cl::Double)) v origin
--      st = showXMDist ( [i|chi2tot ->#{sum cl}, r ->|]) v origin
      sh = showCl (st ++ ", chi2s ->") cl ++ ", Mass ->"
  showMMeas sh $ invMass (map q2p ql)

test :: [String] -> IO ()
test arg =
  case arg of
    ["1"] -> do
          VHMeas v hl <- hSlurp cmsFile
          mapM_ showHelix  hl
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
          showProng $ fit' v hl

-- CMS test file
    ["4"] -> do
          VHMeas v hl <- hSlurp cmsFile
          mapM_ showHelix  hl
          mapM_ showMomentum hl
          doFitTest v hl [0..]

    [fn] -> do
          --mapM_ showMomentum hl
          let
              listMinus1 :: Int -> Int -> [Int]
              listMinus1 n i = filter (/= i) [0..n]
              fitMinus1 :: VHMeas HMeas -> Int -> Prong
              fitMinus1 (VHMeas v hl) = fit' v . hFilter hl . listMinus1 (length hl)

          VHMeas v hl <- hSlurp fn
          mapM_ showMomentum hl
          let nh = length hl - 1
          putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
          mapM_ ( showProng . fitMinus1 (VHMeas v hl) ) [0..nh]

    _ -> exit

doFitTest :: XMeas -> [HMeas] -> [Int] -> IO ()
doFitTest v hl l5 = do
  let showLen xs = show $ length xs
  let showQChi2 (qm, chi2, i) = showQMeas (printf "q%d chi2 ->%8.1f " (i::Int) (chi2::Double) ++ "pt,pz,fi,E ->") qm

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
  putStrLn $ showXMDist (showXMeas "final vertex at" vf ++ ", r =") vf origin


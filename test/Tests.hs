-- file: test.hs
--
module Main ( main ) where

import System.Environment
import System.Exit
import Text.Printf
-- :set -XQuasiQuotes
-- import Data.String.Interpolate

import Input ( hSlurp, dataFiles, hSlurpAll )
import Types (  XMeas (..), HMeas (..), Prong (..), VHMeas (..)
              , DMeas (..), Pos (..), origin
              , invMass, h2p, h2q, q2p
             )
import Fit ( fit, fitw )

import Random ( doRandom )

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

showMomentum :: HMeas -> IO ()
showMomentum h = putStrLn $ "pt,pz,fi,E ->" ++ (show . h2q) h
showHelix :: HMeas -> IO ()
showHelix h = putStrLn $ "Helix ->" ++ (show h)
showProng :: Prong -> IO ()
showProng (Prong _ v ql cl) = do
  let
      showCl :: String -> [Double] -> String
      showCl = foldl (\s x -> s++printf "%8.1g" (x::Double))
      sc = (printf "chi2tot ->%8.1f" (sum cl::Double))
      sd = ", r ->"++ (show $ distance v origin)
      scl = showCl ", chi2s ->" cl
      sm = ", Mass ->" ++ show (invMass (map q2p ql))
  putStrLn $ sc ++ sd ++ scl ++ sm

test :: [String] -> IO ()
test arg =
  case arg of
    ["1"] -> do
          VHMeas v hl <- hSlurp thisFile
          mapM_ showHelix  hl
          mapM_ showMomentum hl
          let l5 = [0,2,3,4,5] -- these are the tracks supposedly from the tau
          doFitTest (VHMeas v hl) l5
          showProng $ fitw (VHMeas v (hFilter hl l5))

    ["2"] -> do
          VHMeas v hl <- hSlurp otherFile
          mapM_ showMomentum hl
          let l5 = [0,1,2,4,5]
          doFitTest (VHMeas v hl) l5
          showProng $ fitw (VHMeas v (hFilter hl l5))

-- slurp in all event data files from ./dat and append helices
    ["3"] -> do
          ps <- dataFiles "dat"
          VHMeas v hl <- hSlurpAll ps
          mapM_ showMomentum hl
          doFitTest (VHMeas v hl) [0..]
          showProng $ fitw (VHMeas v hl)

-- CMS test file
    ["c"] -> do
          VHMeas v hl <- hSlurp cmsFile
--          mapM_ showHelix  hl
--          mapM_ showMomentum hl
--          doFitTest (VHMeas v hl) [0..]
          showProng $ fitw (VHMeas v hl)

    ["r"] -> do
          VHMeas v hl <- hSlurp thisFile
          doRandom 1000 (VHMeas v (hFilter hl [0,2,3,4,5]))

    [fn] -> do
          let
              fitMinus1 :: VHMeas -> Int -> Prong
              fitMinus1 (VHMeas v hl) i = fitw (VHMeas v hl') where
                hl' = hFilter hl (filter (/= i) [0..(length hl)])

          VHMeas v hl <- hSlurp fn
          mapM_ showMomentum hl
          let nh = length hl - 1
          putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
          mapM_ ( showProng . fitMinus1 (VHMeas v hl) ) [0..nh]

    _ -> exit

doFitTest :: VHMeas -> [Int] -> IO ()
doFitTest (VHMeas v hl) l5 = do
  let showLen xs = show $ length xs
  let showQChi2 (qm, chi2, i) = putStrLn $ (printf "q%d chi2 ->%8.1f " (i::Int) (chi2::Double) ++ "pt,pz,fi,E ->") ++ show qm

  putStrLn $ "initial vertex position ->" ++ show v

  let pl              = map h2p hl
  putStrLn $ ("Inv Mass " ++ showLen pl ++ " helix") ++ show (invMass pl)
  let pl5             = map h2p $ hFilter hl l5
  putStrLn $ ("Inv Mass " ++ showLen pl5 ++ " helix") ++ show (invMass pl5)

  putStrLn             "Fitting Vertex --------------------"
  let Prong _ vf ql cl = fit (VHMeas v hl)
  putStrLn $ "Fitted vertex ->" ++ show vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  putStrLn $ "Inv Mass " ++ showLen ql ++ " fit" ++ show (invMass $map q2p ql)
  let pl5              = map q2p $ hFilter ql l5
  putStrLn $ "Inv Mass " ++ showLen pl5 ++ " fit" ++ show (invMass pl5)

  putStrLn             "Refitting Vertex-----------------"
  let Prong _n vf ql cl = fit (VHMeas v (hFilter hl l5))
  putStrLn $ "Refitted vertex ->" ++ show vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  putStrLn $ "Inv Mass " ++ showLen ql ++ " refit" ++ show (invMass $ map q2p ql)
  putStrLn $ ("final vertex at" ++ show vf ++ ", r =") ++ (show $ distance vf origin)


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

-- filter list of objects given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
iflt :: ( Eq a, Enum a, Num a ) => [a] -> [b] -> [b]
iflt rng hl =
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

irem :: (Eq a, Enum a, Num a) => a -> [b] -> [b]
irem indx hl = [ h | (h,i) <- zip hl [0..], i /= indx ]

hFilter :: [Int] -> VHMeas -> VHMeas
hFilter is (VHMeas v hl) = VHMeas v (iflt is hl)

hRemove :: Int -> VHMeas -> VHMeas
hRemove indx (VHMeas v hl) = VHMeas v (irem indx hl)

thisFile  = "dat/tr05129e001412.dat"
otherFile = "dat/tr05158e004656.dat"
thirdFile = "dat/tr05166e001984.dat"
mc0File   = "dat/tr00101e007076.dat"
mc1File   = "dat/tr00101e008340.dat"
mc2File   = "dat/tr00101e012040.dat"
cmsFile   = "dat/tav-0.dat"

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
          vm <- hSlurp thisFile
          mapM_ showHelix $ helices vm
          mapM_ showMomentum $ helices vm
          let l5 = [0,2,3,4,5] -- these are the tracks supposedly from the tau
          doFitTest vm l5
          showProng . fitw . hFilter l5 $ vm

    ["2"] -> do
          vm <- hSlurp otherFile
          mapM_ showMomentum $ helices vm
          let l5 = [0,1,2,4,5]
          doFitTest vm l5
          showProng . fitw . hFilter l5 $ vm

-- slurp in all event data files from ./dat and append helices
    ["3"] -> do
          fs <- dataFiles "dat"
          mapM_ xxx $ drop 1 fs where
            xxx f = do
              vm <- hSlurp $ f
              putStrLn $ printf "File %s" f
              mapM_ showMomentum $ helices vm
              showProng $ fitw vm
              let nh = length (helices vm) - 1
              putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
              mapM_ (\indx -> showProng . fitw . hRemove indx $ vm) [0..nh]

-- slurp in all event data files from ./dat and append helices
    ["4"] -> do
          ps <- dataFiles "dat"
          vm <- hSlurpAll ps
          doFitTest vm [0..]
          showProng $ fitw vm

-- CMS test file
    ["c"] -> do
          vm <- hSlurp cmsFile
--          mapM_ showHelix $ helices vm
          mapM_ showMomentum $ helices vm
          doFitTest vm [0..]
--          showProng $ fitw vm

    ["r"] -> do
--      vm <- hSlurp thisFile
        hSlurp thisFile >>= doRandom 1000 . hFilter [0,2,3,4,5]

    [fn] -> do
          vm <- hSlurp fn
          mapM_ showMomentum $ helices vm
          doFitTest vm [0..]
          let nh = length (helices vm) - 1
          putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
          --mapM_ (\indx -> showProng . fitw . hRemove indx $ vm ) [0..nh]

    _ -> exit

doFitTest :: VHMeas -> [Int] -> IO ()
doFitTest vm l5 = do
  let showLen xs = show $ length xs
  let showQChi2 (qm, chi2, i) = putStrLn $ (printf "q%d chi2 ->%8.1f " (i::Int) (chi2::Double) ++ "pt,pz,fi,E ->") ++ show qm

  putStrLn $ "initial vertex position ->" ++ show (vertex vm)

  let pl              = map h2p $ helices vm
  putStrLn $ ("Inv Mass " ++ showLen pl ++ " helix") ++ show (invMass pl)
  let pl5             = map h2p $ iflt l5 (helices vm)
  putStrLn $ ("Inv Mass " ++ showLen pl5 ++ " helix") ++ show (invMass pl5)

  putStrLn             "Fitting Vertex --------------------"
  let Prong _ vf ql cl = fit vm
  putStrLn $ "Fitted vertex ->" ++ show vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  putStrLn $ "Inv Mass " ++ showLen ql ++ " fit" ++ show (invMass $map q2p ql)
  let pl5              = map q2p $ iflt l5 ql
  putStrLn $ "Inv Mass " ++ showLen pl5 ++ " fit" ++ show (invMass pl5)

  putStrLn             "Refitting Vertex-----------------"
  let Prong _n vf ql cl = fit . hFilter l5 $ vm
  putStrLn $ "Refitted vertex ->" ++ show vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  putStrLn $ "Inv Mass " ++ showLen ql ++ " refit" ++ show (invMass $ map q2p ql)
  putStrLn $ ("final vertex at" ++ show vf ++ ", r =") ++ (show $ distance vf origin)


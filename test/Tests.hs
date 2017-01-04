-- file: test.hs
--
module Main ( main ) where

import System.Environment
import System.Exit
import Text.Printf
import Control.Monad ( liftM )
-- :set -XQuasiQuotes
-- import Data.String.Interpolate

import Input ( hSlurp, dataFiles, hSlurpAll )
import Types (  HMeas (..), Prong (..), VHMeas (..)
              , Pos (..)
  , vBlowup, hFilter, hRemove
              , invMass, h2p, h2q, q2p
             )
import Fit ( fit, fitw, ksm )

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


data DataFileNames = DataFileNames {
    thisFile  :: String
  , otherFile :: String
  , thirdFile :: String
  , badFile   :: String
  , mc0File   :: String
  , mc1File   :: String
  , mc2File   :: String
  , cmsFile   :: String
                                   }

df :: DataFileNames
df = DataFileNames {
    thisFile  = "dat/tr05129e001412.dat"
  , otherFile = "dat/tr05158e004656.dat"
  , thirdFile = "dat/tr05166e001984.dat"
  , badFile   = "dat/tr07849e007984.dat"
  , mc0File   = "dat/tr00101e007076.dat"
  , mc1File   = "dat/tr00101e008340.dat"
  , mc2File   = "dat/tr00101e012040.dat"
  , cmsFile   = "dat/tav-0.dat"
                                   }

showMomentum :: HMeas -> IO ()
showMomentum h = putStrLn $ "pt,pz,fi,E ->" ++ (show . h2q) h
showHelix :: HMeas -> IO ()
showHelix h = putStrLn $ "Helix ->" ++ (show h)
showProng :: Prong -> IO Prong
showProng (Prong n v ql cl m) = do
  let
      showCl :: String -> [Double] -> String
      showCl = foldl (\s x -> s++printf "%8.1g" (x::Double))
      sc = (printf "chi2tot ->%8.1f, ndof %d" (sum cl::Double)) (n*2::Int)
      sd = ", r ->"++ (show $ distance v mempty)
      scl = showCl ", chi2s ->" cl
      sm = ", Mass ->" ++ show (invMass (map q2p ql))
  putStrLn $ sc ++ sd ++ scl ++ sm
  return (Prong n v ql cl m)

test :: [String] -> IO ()
test arg =
  case arg of
    ["1"] -> do
          vm <- liftM (vBlowup 10000.0) (hSlurp . thisFile $ df)

          mapM_ showHelix $ helices vm
          mapM_ showMomentum $ helices vm
          let l5 = [0,2,3,4,5] -- these are the tracks supposedly from the tau
          doFitTest vm l5
          _ <- showProng . fitw . hFilter l5 $ vm
          return ()

    ["2"] -> do
          vm <- liftM (vBlowup 10000.0) (hSlurp . otherFile $ df)
          mapM_ showMomentum $ helices vm
          let l5 = [0,1,2,4,5]
          doFitTest vm l5
          _ <- showProng . fitw . hFilter l5 $ vm
          return ()

    ["3"] -> do -- this file has an additional track that makes the fit bad
          vm <- liftM (vBlowup 10000.0) (hSlurp . badFile $ df)
          mapM_ showHelix $ helices vm
          mapM_ showMomentum $ helices vm
          let l5 = [0,2,3,4,5] -- these are the tracks supposedly from the tau
          doFitTest vm l5
          pr <- showProng . fitw . hFilter l5 $ vm
          let vf = fitVertex $ pr
              h = head . helices . hFilter [6] $ vm
              Just (_, chi2, _) =  ksm vf h
          putStrLn $ printf "chi2 of track 6 w/r to fit vertex is %8.1f" (chi2::Double)
          return ()

-- do tests for each data file
    ["4"] -> do
          fs <- dataFiles "dat"
          mapM_ xxx $ drop 4 fs where
            xxx f = do
              vm <- liftM (vBlowup 10000.0) (hSlurp $ f)
              putStrLn $ printf "File %s" f
              mapM_ showMomentum $ helices vm
              print $ length $ helices vm
              _ <- showProng $ fitw vm
              let nh = length (helices vm) - 1
              putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
              mapM_ (\indx -> showProng . fitw . hRemove indx $ vm) [0..nh]

-- slurp in all event data files from ./dat and append helices
    ["5"] -> do
          ps <- dataFiles "dat"
          vm <- liftM (vBlowup 10000.0) (hSlurpAll ps)
          doFitTest vm [0..]
          _ <- showProng $ fitw vm
          return ()

-- CMS test file
    ["c"] -> do
          vm <- liftM (vBlowup 10000.0) (hSlurp . cmsFile $ df)
--        mapM_ showHelix $ helices vm
          mapM_ showMomentum $ helices vm
          doFitTest vm [0..]
--          showProng $ fitw vm
          return ()

    ["r"] -> do
--          vm <- liftM (vBlowup 10000.0) (hSlurp thisFile)
          (hSlurp . thisFile) df >>= doRandom 1000 . hFilter [0,2,3,4,5] . vBlowup 10000.0

    [fn] -> do
          vm <- liftM (vBlowup 10000.0) (hSlurp fn)
          mapM_ showMomentum $ helices vm
          doFitTest vm [0..]
          let nh = length (helices vm) - 1
          putStrLn $ printf "Inv Mass %d in %d refit, all combinations" (nh::Int) ((nh+1)::Int)
          mapM_ (\indx -> showProng . fitw . hRemove indx $ vm ) [0..nh]

    _ -> exit

doFitTest :: VHMeas -> [Int] -> IO ()
doFitTest vm l5 = do
  let showLen xs = show $ length xs
  let showQChi2 (qm, chi2, i) = putStrLn $ (printf "q%d chi2 ->%8.1f " (i::Int) (chi2::Double) ++ "pt,pz,fi,E ->") ++ show qm

  putStrLn $ "initial vertex position -> " ++ show (vertex vm)

  let pl              = map h2p $ helices vm
  putStrLn $ ("Inv Mass " ++ showLen pl ++ " helix") ++ show (invMass pl)
  let pl5             = map h2p . helices . hFilter l5 $ vm
  putStrLn $ ("Inv Mass " ++ showLen pl5 ++ " helix") ++ show (invMass pl5)

  putStrLn             "Fitting Vertex --------------------"
  let Prong _ vf ql cl _ = fit vm
  putStrLn $ "Fitted vertex -> " ++ show vf
  mapM_ showQChi2 $ zip3 ql cl [0..]
  putStrLn $ "Inv Mass " ++ showLen ql ++ " fit" ++ show (invMass $map q2p ql)

  let m5 = invMass . map q2p . iflt l5 $ ql
      iflt rng hl = [h | (h, i) <- zip hl [0..], i `elem` rng ]
  putStrLn $ "Inv Mass " ++ showLen l5 ++ " fit" ++ show m5

  putStrLn             "Refitting Vertex-----------------"
  let prf = fit . hFilter l5 $ vm
  putStrLn $ "Refitted vertex -> " ++ show (fitVertex prf)
  mapM_ showQChi2 $ zip3 (fitMomenta prf) (fitChi2s prf) [0..]
  putStrLn $ "Inv Mass " ++ (show $ nProng prf) ++ " refit" ++ show (invMass $ map q2p (fitMomenta prf))
  putStrLn $ "Final vertex -> " ++ show (fitVertex prf)


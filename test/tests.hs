-- file: test.hs
--
module Main ( main) where

import Data.Monoid

import Input ( hSlurp, dataFiles )
import Types ( M, V, XMeas (..), HMeas (..), PMeas (..), Prong (..), VHMeas (..)
             , showXMeas, showPMeas, showMMeas )
import Coeff ( w2pt, h2p4, q2p4, invMass )
import Matrix ( inv, sw, fromList, fromList2 )
import Fit ( fit )

hFilter :: ( Eq a, Enum a, Num a ) => [b] -> [a] -> [b]
hFilter hl rng =
-- filter list of tracks of helices etc given list of indices in [a]
-- return list with only those b that have  indices that  are in rng [a]
  [h | (h, i) <- zip hl [0..], i `elem` rng ]

thisFile = "dat/tr05129e001412.dat"
otherFile = "dat/tr05158e004656.dat"
thirdFile = "dat/tr05166e001984.dat"

hSlurpAll :: [String] -> IO (VHMeas HMeas)
hSlurpAll (x:[]) = hSlurp x
hSlurpAll (x:xs) = do
  v0 <- hSlurp x
  vs <- hSlurpAll xs
  return $ v0 <> vs

main :: IO ()
main = do
  VHMeas v hl <- hSlurp thisFile
  showXMeas "initial vertex position ->" v
  mapM_ (showPMeas "px,py,pz,E ->" . h2p4) hl

  let l5 = [0,2,3,4,5]
  doFitTest v hl l5

  VHMeas v hl <- hSlurp otherFile
  mapM_ (showPMeas "px,py,pz,E ->" . h2p4) hl

  let l5 = [0,1,2,4,5]
  doFitTest v hl l5

-- now slurp in all event data files from ./dat and append helices
  ps <- dataFiles "dat"
  VHMeas v hl <- hSlurpAll ps

  mapM_ (showPMeas "px,py,pz,E ->" . h2p4) hl
  doFitTest v hl [0..]
  showXMeas "ok, let's check it" v

doFitTest :: XMeas -> [HMeas] -> [Int] -> IO ()
doFitTest v hl l5 = do
  let showLen xs = show $ length xs

  let pl              = map h2p4 hl
  showMMeas ("Inv Mass " ++ showLen pl ++ " helix") $ invMass pl
  let pl5             = map h2p4 $ hFilter hl l5
  showMMeas ("Inv Mass " ++ showLen pl5 ++ " helix") $ invMass pl5

  putStrLn            "Fitting Vertex --------------------"
  let Prong n vf ql cl = fit v hl
  let pl               = map q2p4 ql
  showMMeas ("Inv Mass " ++ showLen pl ++ " fit") $ invMass pl
  let pl5              = map q2p4 $ hFilter ql l5
  showMMeas ("Inv Mass " ++ showLen pl5 ++ " fit") $ invMass pl5

  putStrLn            "Re-fitting Vertex-----------------"
  let Prong n vf ql cl = fit v $ hFilter hl l5
  let pl               = map q2p4 ql
  showMMeas ("Inv Mass " ++ showLen pl ++ " refit")  $ invMass pl



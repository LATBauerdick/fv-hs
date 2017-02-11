{-# LANGUAGE NoImplicitPrelude, OverloadedStrings #-}

module Cluster ( doCluster ) where

import Types (  VHMeas (..), HMeas (..)
  , XMeas (..), rVertex, zVertex, d0Helix, z0Helix, ptHelix, pzHelix, h2p )
import Fit ( kAdd, ksm )
import Matrix ( sw, scalar, inv )

import Protolude
--import Data.Text
import Text.Printf ( printf )
-- import Control.Monad ( when )
-- import Data.Maybe ( mapMaybe )
--import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

-- import Debug.Trace ( trace )
debug :: a -> [Char] -> a
debug = flip trace


doCluster :: VHMeas -> IO ()
doCluster vm = do
  let v0 = vertex vm -- beamspot
      zs = mapMaybe (zFit v0) $ helices vm
  putStrLn $ "beam spot -> " ++ show v0
  -- let hist = histogram binSturges zs
  let hist = histogramNumBins 200 zs
  _ <- plot "cluster-z.png" hist
  return ()

zFit :: XMeas -> HMeas -> Maybe Double
zFit v h = z where
-- first check chi2 of this helix w/r to vertex position v
  v2 = kAdd v h
  zv2 = zVertex v2
  --  Just (_, chi2', _) =  ksm v h
  XMeas x0 cx0 = v
  XMeas x1 cx1 = v2
  chi2  = scalar . sw (x1-x0) $ inv cx0
  chi2' = scalar . sw (x1-x0) $ inv cx1
  z = if (chi2 >= 0.0) && (chi2 < 1000.0) && (abs zv2) < 10.0
        then Just (zVertex v2)
            `debug` ("vertex track candidate\n" ++ debugVH chi2 chi2' v2 h)
        else Nothing
            `debug` ("failed track candidate\n" ++ debugVH chi2 chi2' v2 h)

  debugVH chi2 chi2' v h =
    printf "-->  chi2=%8.1f (%8.0f)\n" chi2 chi2'
    ++ printf "--> Helix pt=%8.1f GeV, pz=%8.1f GeV, d0=%8.2f cm, z0=%8.2f cm\n"
      (ptHelix h) (pzHelix h) (d0Helix h) (z0Helix h)
    ++ printf "--> Vertex r=%8.2f cm, z=%8.2f cm"
      (rVertex v2) (zVertex v2)


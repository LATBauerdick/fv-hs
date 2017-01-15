module Cluster ( doCluster ) where

import Types (  VHMeas (..), vBlowup, Prong (..), HMeas (..)
  , XMeas (..), zVertex )
import Fit ( fitw, kAdd, ksm )

import Prelude
import Text.Printf ( printf )
import Control.Monad ( when )
import Data.Maybe ( mapMaybe )
import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace


doCluster :: VHMeas -> IO ()
doCluster vm = do
  let v0 = vertex vm -- beamspot
      zs = mapMaybe (ff v0) $ helices vm
  putStrLn $ "beam spot -> " ++ show v0
  -- let hist = histogram binSturges zs
  let hist = histogramNumBins 200 zs
  _ <- plot "cluster-z.png" hist
  return ()

ff :: XMeas -> HMeas -> Maybe Double
ff v h = z where
  v2 = kAdd v h
  zv2 = zVertex v2
  Just (_, chi2, _) =  ksm v h
  z = if (chi2 >= 0.0) && (chi2 < 1000.0) && (abs zv2) < 10.0
        then Just (zVertex v2)
            `debug` ("vertex track candidate " ++ (printf "chi2 %6.0f, v at " chi2) ++ (show v2))
        else Nothing

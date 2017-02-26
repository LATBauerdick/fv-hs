{-# LANGUAGE NoImplicitPrelude, OverloadedStrings #-}

module Cluster ( doCluster ) where

import Types (  VHMeas (..), HMeas (..), Prong (..)
  , XMeas (..), XFit (..)
  , rVertex, chi2Vertex, zVertex, d0Helix, z0Helix, ptHelix, pzHelix, h2p )
import Fit ( kAdd, kAddF, ksm )
import Matrix ( sw, scalar, inv )
import Coeff ( qv2h, gamma )

--import Protolude
--import Data.Text

import Prelude
import Data.Maybe ( mapMaybe )
import Data.List ( sortOn )

import Text.Printf ( printf )
-- import Control.Monad ( when )
-- import Data.Maybe ( mapMaybe )
--import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

import Debug.Trace ( trace )
debug :: a -> [Char] -> a
debug = flip trace


doCluster :: VHMeas -> IO ()
doCluster vm = do
  let v0 = vertex vm -- beamspot
      zs = mapMaybe (zFit (xFit v0)) $ helices vm
  putStrLn $ "beam spot -> " ++ show v0
  -- let hist = histogram binSturges zs
  let hist = histogramNumBins 40 zs
  _ <- plot "cluster-z.png" hist
  putStrLn "---------------------------------------------------"
  putStrLn "---------------------------------------------------"
  putStrLn "---------------------------------------------------"
  putStrLn "---------------------------------------------------"
  print . vTree . cleanup $ vm
  return ()

xFit :: XMeas -> XFit
xFit (XMeas v vv) = XFit v vv 1e6
xMeas :: XFit -> XMeas
xMeas (XFit v vv _) = XMeas v vv
cleanup :: VHMeas -> VHMeas
-- remove vm helices that are incompatible with vm vertex
cleanup (VHMeas v hl) = (VHMeas v hl') where
  hl' = sortOn z0Helix . mapMaybe (fff v) $ hl
  fff :: XMeas -> HMeas -> Maybe HMeas
  fff v h = mh where
-- first check chi2 of this helix w/r to vertex position v
    vf  = kAddF (xFit v) h
    zvf = zVertex vf
    chi2f = chi2Vertex vf
    XMeas x0 cx0 = v
    XFit x1 cx1 _ = vf
    chi2  = scalar . sw (x1-x0) $ inv cx0
    chi2' = scalar . sw (x1-x0) $ inv cx1
    good = (chi2f < 1000.0) && (abs zvf) < 10.0
      `debug` (printf "---> chi2=%8.1f chi2=%8.1f chi2=%8.1f gamma=%9.3f\n" chi2f chi2 chi2' (180.0/pi*(Coeff.gamma h (xMeas vf))))
    mh = if good
          then Just h
            `debug` ("h does belong to beam spot \n" ++ debugVH chi2 chi2' (xMeas vf) h)
          else Nothing
            `debug` ("h does not belong to beam spot\n" ++ debugVH chi2 chi2' (xMeas vf) h)

debugVH :: Double -> Double -> XMeas -> HMeas -> String
debugVH chi2 chi2' v h =
    printf "-->  chi2=%8.1f (%8.0f)\n" chi2 chi2'
    ++ printf "--> Helix pt=%8.1f GeV, pz=%8.1f GeV, d0=%8.2f cm, z0=%8.2f cm\n"
      (ptHelix h) (pzHelix h) (d0Helix h) (z0Helix h)
    ++ printf "--> Vertex r=%8.2f cm, z=%8.2f cm"
      (rVertex (xFit v)) (zVertex (xFit v))


data HTree a = Empty | Node a (HTree a) deriving (Show)

vTree :: VHMeas -> HTree Prong
vTree vm = Node p vRight where
  (p,vmr) = cluster vm
  vRight = case vmr of
             Nothing -> Empty
             Just vm' -> vTree vm'


cluster :: VHMeas -> (Prong, Maybe VHMeas)
cluster (VHMeas v hl) = ( p, r ) where
  p = kSmooth (VHMeas v hl) . foldl kAdd v $ hl
  r = Nothing

kSmooth :: VHMeas -> XMeas -> Prong
--kSmooth vm v | trace ("kSmooth " ++ (show . length . view helicesLens $ vm) ++ ", vertex at " ++ (show v) ) False = undefined
kSmooth (VHMeas v0 hl) v = pr' where
  (ql, chi2l, hl') = unzip3 $ mapMaybe (ksm v) hl
  (n, n') = (length hl, length ql)
  n'' = if n == n' then n else n' `debug` "kSmooth killed helices"
  pr' = Prong { fitVertex = v, fitMomenta = ql, fitChi2s = chi2l, nProng = n'', measurements = VHMeas v0 hl' }

{-
fitCluster :: VHMeas -> IO ()
fitCluster vm = do
  let v0 = vertex vm
      hl = helices vm
  foldl addOne v0 hl

addOne :: XMeas -> HMeas -> Maybe (Qmeas, Chi2, HMeas)
addOne v h t = do -- add only if chi2 at temperature t is reasonable
  let v' = kAddF v h
      XMeas x0 cx0 = v
      XMeas x1 cx1 = v'
      chi2  = scalar . sw (x1-x0) $ inv cx0
      w = wght t chi2
      v'' = if (w > 0.9)
               then v'
               else v

  Just (q, chi2, _) <- ksm v' h
  if (chi2 < cut) then return Just (q, chi2, h)
                  else return Nothing
-}


zFit :: XFit -> HMeas -> Maybe Double
zFit v h = z where
-- first check chi2 of this helix w/r to vertex position v
  v2 = kAddF v h
  zv2 = zVertex v2
  --  Just (_, chi2', _) =  ksm v h
  XMeas x0 cx0 = xMeas v
  XMeas x1 cx1 = xMeas v2
  chi2  = scalar . sw (x1-x0) $ inv cx0
  chi2' = scalar . sw (x1-x0) $ inv cx1
  z = if (chi2 >= 0.0) && (chi2 < 1000.0) && (abs zv2) < 10.0
        then Just (zVertex v2)
            `debug` ("vertex track candidate\n" ++ debugVH chi2 chi2' (xMeas v2) h)
        else Nothing
            `debug` ("failed track candidate\n" ++ debugVH chi2 chi2' (xMeas v2) h)


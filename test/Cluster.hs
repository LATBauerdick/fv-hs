{-# LANGUAGE NoImplicitPrelude, OverloadedStrings #-}

module Cluster ( doCluster ) where

import Types (  VHMeas (..), HMeas (..), QMeas (..), Prong (..)
  , XMeas (..), XFit (..)
  , rVertex, chi2Vertex, zVertex, d0Helix, z0Helix, ptHelix, pzHelix, h2p )
import Fit ( kAdd, kAddF, ksm )
import Matrix ( sw, scalar, inv )
import Coeff ( qv2h, gamma )

--import Protolude
--import Data.Text

import Prelude
import Data.Maybe ( mapMaybe, isNothing, isJust )
import Data.List ( sortOn )
import qualified Math.Gamma ( q )

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
  putStrLn $ "beam spot -> " ++ show v0

  let histz = histogramNumBins 10 $ zs vm
  _ <- plot "cluster-z.png" histz

  let histp = histogramNumBins 11 $ 1.0 : 0.0 : probs vm
  _ <- plot "pd.png" histp

  let Node p0 ht = vTree $ (cleanup vm)
  print $ fitChi2s p0
  print $ fitVertex p0
  print $ nProng p0
  case ht of
    Empty     -> putStrLn "Empty"
    Node p1 _ -> print $ fitVertex p1
  return ()

xFit :: XMeas -> XFit
xFit (XMeas v vv) = XFit v vv 1e6
xMeas :: XFit -> XMeas
xMeas (XFit v vv _) = XMeas v vv

probs :: VHMeas -> [Double]
probs (VHMeas v hl) = filter (\x -> x>0.01) $ map (Math.Gamma.q 1.0 . chi2Vertex . kAddF (xFit v)) hl
zs :: VHMeas -> [Double]
zs (VHMeas v hl) = filter (\x -> (abs x)<10.0) $ map (zVertex . kAddF (xFit v)) hl

cleanup :: VHMeas -> VHMeas
-- remove vm helices that are incompatible with vm vertex
cleanup (VHMeas v hl) = (VHMeas v hl') where
  hl' = sortOn z0Helix . mapMaybe (fff v) $ hl
  fff :: XMeas -> HMeas -> Maybe HMeas
  fff v h = mh where
-- check chi2 of this helix w/r to vertex position v
    vf  = kAddF (xFit v) h
    zvf = zVertex vf
    chi2f = chi2Vertex vf
    prob = Math.Gamma.q 1.0 (chi2f/2.0) -- chi2 distribution with NDOF=2
    good = (prob > 0.01) && (abs zvf) < 10.0
    mh = if good
            -- `debug` (printf "--> chi2=%8.1f prob=%8.4f z=%9.3f gamma=%9.3f" chi2f prob zvf (180.0/pi*(Coeff.gamma h (xMeas vf))))
          then Just h
          else Nothing

data HTree a = Empty | Node a (HTree a) deriving (Show)

vTree :: VHMeas -> HTree Prong
vTree vm = Node p vRight where
  (p,vmr) = cluster vm
  vRight = case vmr of
             Nothing -> Empty
             Just vm' -> vTree vm'


cluster :: VHMeas -> (Prong, Maybe VHMeas)
cluster (VHMeas v hl) = ( p, r ) where
  v1 = foldl kAdd v hl
  p = ks (VHMeas v hl) . foldl kAdd v $ hl
  r = Nothing

--ks vm v | trace ("kSmooth " ++ (show . length . view helicesLens $ vm) ++ ", vertex at " ++ (show v) ) False = undefined
ks (VHMeas v0 hl) v = pr' where
  ll :: [(Maybe (QMeas, Double, HMeas), HMeas)]
  ll = zip (map (ksm v) hl) hl
  hlnothing = [ h | (k, h) <- ll, isNothing k]
  hljust    = [ h | (k, h) <- ll, isJust k]
  (ql, chi2l, hl') = unzip3 $ mapMaybe (ksm v) hl
  (n, n') = (length hl, length ql)
  n'' = if n == n' then n else n' `debug` "kSmooth killed helices"
  pr' = Prong { fitVertex = v, fitMomenta = ql, fitChi2s = chi2l, nProng = n'', measurements = VHMeas v0 hl' }

-- gamma function P(a, x) from https://wiki.haskell.org/Gamma_and_Beta_function
-- approximation is taken from [Lanczos, C. 1964 SIAM Journal on Numerical Analysis,
-- ser. B, vol. 1, pp. 86-96]

-- standard normal CDF https://www.johndcook.com/blog/haskell-erf/
erf :: Double -> Double
erf x = sign*y
    where
        a1 =  0.254829592
        a2 = -0.284496736
        a3 =  1.421413741
        a4 = -1.453152027
        a5 =  1.061405429
        p  =  0.3275911

        -- Abramowitz and Stegun formula 7.1.26
        sign = if x > 0
                   then  1
                   else -1
        t  =  1.0/(1.0 + p* abs x)
        y  =  1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)

test_erf :: Bool
test_erf = maximum [ abs(erf x - y)  | (x, y) <- zip xs ys ] < epsilon
    where
        epsilon = 1.5e-7 -- accuracy promised by A&S
        xs = [-3, -1, 0.0, 0.5, 2.1 ]
        ys = [-0.999977909503,
              -0.842700792950,
               0.0,
               0.520499877813,
               0.997020533344]

-- standard normal CDF https://www.johndcook.com/blog/haskell-phi/
phi :: Double -> Double
phi x = y 
    where
        a1 =  0.254829592
        a2 = -0.284496736
        a3 =  1.421413741
        a4 = -1.453152027
        a5 =  1.061405429
        p  =  0.3275911

        -- Abramowitz and Stegun formula 7.1.26
        sign = if x > 0
                   then  1
                   else -1
        t = 1.0/(1.0 + p * abs x / sqrt 2.0)
        e = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x/2.0)
        y = 0.5*(sign*e + 1.0)

test_phi :: Bool
test_phi = maximum [ abs(phi x - y) | (x, y) <- zip xs ys ] < epsilon
    where
        epsilon = 1.5e-7 -- accuracy promised by A&S
        xs = [-3, -1, 0.0, 0.5, 2.1 ]
        ys = [0.00134989803163, 
              0.158655253931, 
              0.5, 
              0.691462461274, 
              0.982135579437]


-- returns the number of samples and the mean,
-- but does not directly return variance, skewness and kurtosis.
-- Instead it returns moments from which these statistics can easily be calculated
-- using the mvks function
moments (n, m1, m2, m3, m4) x = (n', m1', m2', m3', m4')
        where
            n' = n + 1
            delta = x - m1
            delta_n = delta / n'
            delta_n2 = delta_n**2
            t = delta*delta_n*n
            m1' = m1 + delta_n
            m4' = m4 + t*delta_n2*(n'*n' - 3*n' + 3) + 6*delta_n2*m2 - 4*delta_n*m3
            m3' = m3 + t*delta_n*(n' - 2) - 3*delta_n*m2
            m2' = m2 + t

mvsk (n, m1, m2, m3, m4) = (m1, m2/(n-1.0), (sqrt n)*m3/m2**1.5, n*m4/m2**2 - 3.0)




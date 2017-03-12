{-# LANGUAGE NoImplicitPrelude, OverloadedStrings #-}

module FVT.Cluster ( doCluster, fsmw ) where

import FV.Types (  VHMeas (..), HMeas (..), QMeas (..), Prong (..)
  , XMeas (..), XFit (..)
  , chi2Vertex, zVertex, z0Helix )
import FV.Fit ( kAdd, kAddF, ksm )
import FV.Matrix ( fromList, toList )
import FV.Types ( X3 )

--import Data.Text

import Prelude
import Data.Maybe ( mapMaybe, isNothing, isJust )
import Data.List ( sortOn, foldl' )
import qualified Math.Gamma ( q )

import Text.Printf ( printf )
-- import Control.Monad ( when )
-- import Data.Maybe ( mapMaybe )
-- import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

doCluster :: VHMeas -> IO ()
doCluster vm = do
  let v0 = vertex vm -- beamspot
  putStrLn $ "beam spot -> " ++ show v0

  let histz = histogramNumBins 20 $ zs vm
  _ <- plot "cluster-z.png" histz

  let histp = histogramNumBins 11 $ 1.0 : 0.0 : probs vm
  _ <- plot "pd.png" histp

  let Node p0 ht = vTree $ cleanup vm
  print $ fitChi2s p0
  print $ fitVertex p0
  print $ nProng p0
  case ht of
    Empty     -> putStrLn "Empty"
    Node p1 _ -> print $ fitVertex p1
  return ()

xFit :: XMeas -> XFit
xFit (XMeas v vv) = XFit v vv 1e6

probs :: VHMeas -> [Double]
probs (VHMeas v hl) = filter (\x -> x>0.01) $ map (Math.Gamma.q 1.0 . chi2Vertex . kAddF (xFit v)) hl
zs :: VHMeas -> [Double]
zs (VHMeas v hl) = filter (\x -> (abs x)<10.0) $ map (zVertex . kAddF (xFit v)) hl

cleanup :: VHMeas -> VHMeas
-- remove vm helices that are incompatible with vm vertex
cleanup (VHMeas v hl) = (VHMeas v hl') where
  hl' = sortOn z0Helix . mapMaybe (filtProb 0.01 v) $ hl

filtProb :: Double -> XMeas -> HMeas -> Maybe HMeas
filtProb cut v h = mh where
-- check chi2 of this helix w/r to vertex position v
  vf  =   kAddF (xFit v) h
  zvf =   zVertex vf
  chi2f = chi2Vertex vf
  prob =  Math.Gamma.q 1.0 (chi2f/2.0) -- chi2 distribution with NDOF=2
  good =  (prob > cut) && (abs zvf) < 10.0
  mh =    if good
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
  zs  = filter (\x -> (abs x)<10.0) . map (zVertex . kAddF (xFit v)) $ hl
  z0 = fsmw (length zs) zs
  XMeas x0 cx0 = v
  [xx0, yx0] = FV.Matrix.toList 2 x0
  v0 = XMeas (FV.Matrix.fromList 3 [xx0, yx0, z0]) cx0 `debug` (printf "--> cluster z0=%9.3f " z0)
  v1        = foldl kAdd v0 hl `debug` (show v0)
  ll        = zip (map (ksm v) hl) hl
  hlnothing = [ h | (k, h) <- ll, isNothing k]
  hljust    = [ h | (k, h) <- ll, isJust k]
  qljust    = [ q | (Just (q, _, _), _) <- ll]
  c2just    = [ c | (Just (_, c, _), _) <- ll]
  hlfilt    = mapMaybe (filtProb 0.00000000000000000001 v1) hl `debug` (printf "--> nothing=%5d just=%5d %5d" (length hlnothing) (length hljust) (length qljust))

  (ql, chi2l, hl') = unzip3 $ mapMaybe (ksm v1) hlfilt
  p = Prong { fitVertex = v1, fitMomenta = ql, fitChi2s = chi2l, nProng = (length ql), measurements = VHMeas v hl' }
--  p = ks (VHMeas v hl) v1
--  r = Nothing

  p000 = Prong { fitVertex = v1, fitMomenta = qljust, fitChi2s = c2just, nProng = (length c2just), measurements = VHMeas v hljust }
  r = case (length hlnothing) of
        0 -> Nothing
        _ -> Just (VHMeas v hlnothing)

--ks vm v | trace ("kSmooth " ++ (show . length . view helicesLens $ vm) ++ ", vertex at " ++ (show v) ) False = undefined

-- -------------------------------------------------------
-- use robust Mode finding to get initial vertex position
-- -------------------------------------------------------
-- Fraction-of Sample Mode with Weight method,
-- see Frühwirth & Waltenberger CMS Note 2007/008
--newtype WeightedPoint = WeightedPoint Double
type WeightedPoint = Double
fsmw :: Int -> [WeightedPoint] -> WeightedPoint
fsmw 1 [x] = x
fsmw 2 [x0, x1] = 0.5 * (x1+x0)
fsmw 3 [x0, x1, x2] = case ( 2*x1-x0-x2 ) of
                       xx | xx < 0 -> (x0+x1)/2.0
                       xx | xx > 0 -> (x1+x2)/2.0
                       xx | xx == 0 -> x1
-- fsmw 4 [x0, x1, x2, x3] = x where
--   w  = weightOfInterval [x0, x1, x2, x3]
--   w0 = weightOfInterval [x0, x1, x2]
--   w1 = weightOfInterval [x1, x2, x3]
--   x  = if w0/w < w1/w then fsmw 3 [x0, x1, x2] else fsmw 3 [x1, x2, x3] `debug` (show w ++ ", " ++ show w0 ++ ", " ++ show w1)
fsmw n xs = h where
  h = if n>6 then hsm n xs -- for large n, fall back to no-weight hsm formula
          else fsmw n' xs' where
            alpha :: Double; alpha = 0.5
            n' = ceiling (fromIntegral n * alpha)
            findMin :: (WeightedPoint, Int) -> (WeightedPoint, Int) -> (WeightedPoint, Int)
            findMin (w0, j0) (w, j) = if w<w0 || w0<0 then (w, j) else (w0, j0)
            (_, j') = foldl' findMin (-1.0, 0) . zip (map ( \j ->  weightOfInterval . take n' . drop j $ xs ) $ [0..n-n']) $ [0..]
            xs' = take n' . drop (j') $ xs -- `debug` ("fsmw--> " ++ show n ++ ", " ++ show n' ++ ", " ++ show j' ++ ", " ++ show xs)

-- HalfSampleMode 
-- see Bickel & Frühwirth, On a Fast, Robust Estimator of the Mode, 2006
-- http://ideas.repec.org/a/eee/csdana/v50y2006i12p3500-3530.html
hsm :: Int -> [WeightedPoint] -> WeightedPoint
hsm n xs = fsmw n' xs' where
  alpha :: Double; alpha = 0.5
  n' = ceiling (fromIntegral n * alpha)
  xns = drop (n'-1) xs
  wmin = (last xs - head xs)
  findMin :: (WeightedPoint, Int) -> (WeightedPoint, Int) -> (WeightedPoint, Int)
  findMin (w0, j0) (w, j) = if w<w0 then (w, j) else (w0, j0)
  calcWeight :: WeightedPoint -> WeightedPoint -> WeightedPoint
  calcWeight = (-)
  (_, j') = foldl' findMin (wmin, 0) . zip (zipWith calcWeight xns xs) $ [0..]
  xs' = take n' . drop (j') $ xs --`debug` ("hsm---> " ++ show n ++ ", " ++ show n' ++ ", " ++ show j' ++ ", " ++ show xs ++ ", " ++ show xns)

wDist :: WeightedPoint -> WeightedPoint -> WeightedPoint
wDist x0 x1 = 1.0 / sqrt (x1 - x0 + dmin) where dmin = 0.001 -- 10 µm
weightOfInterval :: [WeightedPoint] -> WeightedPoint
weightOfInterval xs = w where --`debug` ("weightOfInterval--> " ++ show xs ++ ", " ++ show ws ++ ", " ++ show w) where
  ws = [ wDist x0 x1 | x0 <- xs , x1 <- xs, x1>x0]
  w  = (last xs - head xs) / sum ws




{-
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
-}


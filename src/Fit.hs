-- file src/Fit.hs
module Fit ( fit, fitw ) where

import Types (  XMeas (..), HMeas (..), QMeas (..), VHMeas (..)
              , Prong (..), Jaco (..), Chi2
              , X3, C33, Q3
             )
import qualified Coeff ( expand, hv2q )
import Matrix ( inv, tr, sw, scalar, scale)
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

wght :: Double -> Chi2 -> Double -- weight function with Temperature t
wght t chi2 = w where
  chi2cut = 9.0
  w = 1.0/(1.0 + exp ((chi2-chi2cut)/2.0/t))

fitw :: VHMeas -> Prong -- fit with annealing function
fitw (VHMeas v0 hl) = pr where
  -- VHMeas v _ = kFilter (VHMeas v0 hl)
  Prong _ _ _ cl = kSmooth . kFilter $ (VHMeas v0 hl)
  ws  = fmap (wght 10.0) cl
  Prong _ _ _ cl' = fitww ws (VHMeas v0 hl)
  ws' = fmap (wght 1.0) cl'
  pr  = fitww ws' (VHMeas v0 hl)

fitww :: [Double] -> VHMeas -> Prong
fitww ws = kSmoothW . (kFilterW ws)

kFilterW :: [Double] -> VHMeas -> VHMeas
kFilterW ws (VHMeas x ps) = VHMeas x' ps' where
  invWght :: HMeas -> Double -> HMeas
  invWght (HMeas h hh w0) wght = HMeas h (scale wght (inv hh)) w0
  ps' = zipWith invWght ps ws
  x' = foldl kAddW x ps'

kAddW :: XMeas -> HMeas -> XMeas
kAddW (XMeas v vv) (HMeas h gg w0) = kAdd' x_km1 p_k x_e q_e 1e6 0 where
  x_km1 = XMeas v (inv vv)
  p_k   = HMeas h gg w0
  x_e   = v
  q_e   = Coeff.hv2q h v

kSmoothW :: VHMeas -> Prong
kSmoothW (VHMeas v ps') = Prong (length ql) v ql chi2l where
  (ql, chi2l) = unzip $ map (ksm' v) ps'

goodEnough :: Double -> Double -> Int -> Bool
goodEnough c0 c i = abs (c - c0) < chi2cut || i > iterMax where
  chi2cut = 0.5
  iterMax = 99 :: Int

-- kalman smooth: calculate 3-mom q at kalman filter vertex v
ksm' :: XMeas -> HMeas -> (QMeas, Chi2)
ksm' (XMeas x cc) (HMeas h gg w0) = (QMeas q dd w0, chi2) where
    Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
    aaT  = tr aa-- (aa ^+)
    bbT  = tr bb
    ww   = inv (sw bb gg)
    p    = h - h0
    uu   = inv cc
    q    = ww * bbT * gg * (p - aa * x)
    ee   = - cc * aaT * gg * bb * ww
    dd   = ww + sw ee uu
    r    = p - aa*x - bb*q
    gb   = gg - sw gg (sw bbT ww)
    uu'  =  uu - sw aa gb
    cc'  = inv uu'
    x'   = cc' * (uu*x - aaT * gb * p)
    dx   = x - x'
    chi2 = scalar (sw dx uu' + sw r gg)

{-
  data Prong = Prong N XMeas [QMeas] [Chi2] ...
  data VHMeas = VHMeas XMeas [HMeas] ...
  instance Monoid VHMeas where ...
-}
fit :: VHMeas -> Prong
fit = kSmooth . kFilter

kFilter :: VHMeas -> VHMeas
kSmooth :: VHMeas -> Prong
kFilter (VHMeas x ps) = VHMeas (foldl kAdd x ps) ps

kAdd :: XMeas -> HMeas -> XMeas
kAdd (XMeas v vv) (HMeas h hh w0) = kAdd' x_km1 p_k x_e q_e 1e6 0 where
  x_km1 = XMeas v (inv vv)
  p_k   = HMeas h (inv hh) w0
  x_e   = v
  q_e   = Coeff.hv2q h v

kAdd' :: XMeas -> HMeas -> X3 -> Q3 -> Double -> Int -> XMeas
kAdd' (XMeas v0 uu0) (HMeas h gg w0) x_e q_e 𝜒2_0 iter = x_k where
  Jaco aa bb h0 = Coeff.expand x_e q_e
  aaT   = tr aa; bbT = tr bb
  ww    = inv (sw bb gg)
  gb    = gg - sw gg (sw bbT ww)
  uu    = uu0 + sw aa gb; cc = inv uu
  m     = h - h0
  v     = cc * (uu0 * v0 + aaT * gb * m)
  dm    = m - aa * v
  q     = ww * bbT * gg * dm
  𝜒2    = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
  x_k   = if goodEnough 𝜒2_0 𝜒2 iter
            then XMeas v cc
            else kAdd' (XMeas v0 uu0) (HMeas h gg w0) v q 𝜒2 (iter+1)

kSmooth (VHMeas v hl) = Prong (length ql) v ql chi2l where
  (ql, chi2l) = unzip $ map (ksm v) hl

-- kalman smooth: calculate 3-mom q and chi2 at kalman filter vertex
ksm :: XMeas -> HMeas -> (QMeas, Chi2)
ksm (XMeas x cc) (HMeas h hh w0) = (QMeas q dd w0, chi2) where
    Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
    aaT  = tr aa
    bbT  = tr bb
    gg   = inv hh
    ww   = inv (sw bb gg)
    p    = h - h0
    uu   = inv cc
    q    = ww * bbT * gg * (p - aa * x)
    ee   = - cc * aaT * gg * bb * ww
    dd   = ww + sw ee uu
    r    = p - aa*x - bb*q
    gb   = gg - sw gg (sw bbT ww)
    uu'  =  uu - sw aa gb
    cc'  = inv uu'
    x'   = cc' * (uu*x - aaT * gb * p)
    dx   = x - x'
    chi2 = scalar (sw dx uu' + sw r gg)


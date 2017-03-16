-- file src/Fit.hs
module FV.Fit ( fit, fitw, ksm, ksm', kAddF, kAdd ) where

import FV.Types (  XMeas (..), HMeas (..), QMeas (..), VHMeas (..), XFit (..)
--  , helicesLens, view, over, set
              , Prong (..), Jaco (..), Chi2
              , X3, Q3
             )

import qualified FV.Coeff as Coeff ( expand, hv2q )
import FV.Matrix ( inv, invMaybe, det, tr, sw, scalar, scale)

import Prelude
import Data.Maybe ( mapMaybe )

import Text.Printf
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

wght :: Double -> Chi2 -> Double -- weight function with Temperature t
wght t chi2 = w where
  chi2cut = 9.0
  w = 1.0/(1.0 + exp ((chi2-chi2cut)/2.0/t))

fit :: VHMeas -> Prong
fit vm = kSmooth vm . kFilter $ vm

fitw :: VHMeas -> Prong -- fit with annealing function
fitw vm = pr where
  ws  = fmap (wght 10.0) $ fitChi2s . kSmooth vm . kFilter $ vm
  -- ws' = fmap (wght  1.0) $ fitChi2s . kSmooth vm . kFilterW ws $ vm
  ws' = fmap (wght  1.0) $ fitChi2s . kSmooth vm . kFilterW (repeat 1.0) $ vm
          `debug` ("fitw: " ++ (foldl (\ z a -> z ++ printf "%8.3f, " (a::Double))) "weights-> " ws)
  -- pr  = kSmooth vm . kFilterW ws' $ vm
  pr  = kSmooth vm . kFilterW (repeat 1.0) $ vm
          `debug` ("fitw: " ++ (foldl (\ z a -> z ++ printf "%8.3f, " (a::Double))) "weights-> " ws')

kFilter :: VHMeas -> XMeas
kFilter (VHMeas v hl) = foldl kAdd v hl

kAdd :: XMeas -> HMeas -> XMeas
kAdd (XMeas v vv) (HMeas h hh w0) = kAdd' x_km1 p_k x_e q_e 1e6 0 where
  x_km1 = XMeas v (inv vv)
  p_k   = HMeas h (inv hh) w0
  x_e   = v
  q_e   = Coeff.hv2q h v

kAddF :: XFit -> HMeas -> XFit
kAddF (XFit v vv _) (HMeas h hh w0) = kAddF' v (inv vv) h (inv hh) v (Coeff.hv2q h v) 1e6 0

kAddF' v0 uu0 h gg x_e q_e 𝜒2_0 iter = x_k where
    Jaco aa bb h0 = Coeff.expand x_e q_e
    aaT   = tr aa; bbT = tr bb
    x_k   = case invMaybe (sw bb gg) of
              Nothing  -> (XFit v0 (inv uu0) 1e6)  `debug` "... in kAdd'"
              Just ww' -> x_k' where 
                ww    = ww'
                gb    = gg - sw gg (sw bbT ww)
                uu    = uu0 + sw aa gb; cc = inv uu
                m     = h - h0
                v     = cc * (uu0 * v0 + aaT * gb * m)
                dm    = m - aa * v
                q     = ww * bbT * gg * dm
                𝜒2    = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
                x_k'  = if goodEnough 𝜒2_0 𝜒2 iter -- `debug` (printf "--> kAddF' chi2 is %9.1f, %9.1f" 𝜒2 (scalar $ sw (v-v0) uu0))
                          then XFit v cc 𝜒2
                          else kAddF' v0 uu0 h gg v q 𝜒2 (iter+1)

kFilterW :: [Double] -> VHMeas -> XMeas
kFilterW ws vm = v' where
  VHMeas v hl = vm
  v' = foldl kAddW v $ zip hl ws

kAddW :: XMeas -> (HMeas, Double) -> XMeas
kAddW (XMeas v vv) (hm, w) = kAdd' x_km1 p_k x_e q_e 1e6 0 where
  x_km1 = XMeas v (inv vv)
  HMeas h hh w0 = hm
  p_k   = HMeas h (scale w (inv hh)) w0
  x_e   = v
  q_e   = Coeff.hv2q h v

goodEnough :: Double -> Double -> Int -> Bool
--goodEnough c0 c i | trace ("."++show i ++ "|" ++ printf "%8.1f" (abs (c-c0)) ++ printf "%8.1f" c) False = undefined
goodEnough c0 c i = abs (c - c0) < chi2cut || i > iterMax where
  chi2cut = 0.5
  iterMax = 99 :: Int

-- add a helix measurement to kalman filter, return updated vertex position
-- if we can't invert, don't update vertex
kAdd' :: XMeas -> HMeas -> X3 -> Q3 -> Double -> Int -> XMeas
kAdd' (XMeas v0 uu0) (HMeas h gg w0) x_e q_e 𝜒2_0 iter = x_k where
  Jaco aa bb h0 = Coeff.expand x_e q_e
  aaT   = tr aa; bbT = tr bb
  x_k   = case invMaybe (sw bb gg) of
            Just ww' -> x_k' where 
              ww    = ww'
              gb    = gg - sw gg (sw bbT ww)
              uu    = uu0 + sw aa gb; cc = inv uu
              m     = h - h0
              v     = cc * (uu0 * v0 + aaT * gb * m)
              dm    = m - aa * v
              q     = ww * bbT * gg * dm
              𝜒2    = scalar $ sw (dm - bb * q) gg + sw (v - v0) uu0
              x_k'  = if goodEnough 𝜒2_0 𝜒2 iter -- `debug` ("--> kAdd' chi2 is " ++ show 𝜒2)
                then XMeas v cc
                else kAdd' (XMeas v0 uu0) (HMeas h gg w0) v q 𝜒2 (iter+1)
            Nothing -> (XMeas v0 (inv uu0))  `debug` "... in kAdd'"

kSmooth :: VHMeas -> XMeas -> Prong
--kSmooth vm v | trace ("kSmooth " ++ (show . length . view helicesLens $ vm) ++ ", vertex at " ++ (show v) ) False = undefined
kSmooth (VHMeas v0 hl) v = pr' where
  (ql, chi2l, hl') = unzip3 $ mapMaybe (ksm v) hl
  (n, n') = (length hl, length ql)
  n'' = if n == n' then n else n' `debug` "kSmooth killed helices"
  pr' = Prong { fitVertex = v, fitMomenta = ql, fitChi2s = chi2l, nProng = n'', measurements = VHMeas v0 hl' }

-- kalman smoother step: calculate 3-mom q and chi2 at kalman filter'ed vertex
-- if we can't invert, return Nothing and this track will not be included
ksm :: XMeas -> HMeas -> Maybe (QMeas, Chi2, HMeas)
ksm (XMeas x cc) hm = do
  let
      HMeas h hh w0 = hm
      Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
      gg   = inv hh
  ww <- invMaybe (sw bb gg)
  let
      p    = h - h0
      uu   = inv cc
      aaT  = tr aa; bbT   = tr bb
      q    = ww * bbT * gg * (p - aa * x)
      ee   = - cc * aaT * gg * bb * ww
      dd   = ww + sw ee uu
      r    = p - aa*x - bb*q
      ch   = scalar $ sw r gg

      gb   = gg - sw gg (sw bbT ww)
      uu'  =  uu - sw aa gb
      duu  = det uu'
      bad  = duu < 0
      cx   = if bad then 1000.0 `debug` ("--> ksm bad" ++ show duu ++ show uu')
                    else cx'' where
                      cc'  = inv uu' -- `debug` ("--> ksm " ++ show uu')
                      x'   = cc' * (uu*x - aaT * gb * p)
                      dx   = x - x'
                      cx'  = scalar $ sw dx uu'
                      cx'' = if cx' < 0 then 2000.0 `debug` ("--> ksm chi2 is " ++ show cx' ++ ", " ++ show ch ++ ", " ++ show ((max cx' 0) + ch))
                                        else cx'
      chi2 = cx + ch
  return (QMeas q dd w0, chi2, hm)

-- kalman smoother step: calculate 3-mom q and chi2 at kalman filter'ed vertex
-- if we can't invert, return Nothing and this track will not be included
ksm' :: XMeas -> Maybe HMeas -> Maybe (QMeas, Chi2)
ksm' _ Nothing = Nothing
ksm' (XMeas x cc) (Just (HMeas h hh w0)) = do
  let
      Jaco aa bb h0 = Coeff.expand x (Coeff.hv2q h x)
      gg   = inv hh
  ww <- invMaybe (sw bb gg)
  let
      p    = h - h0
      uu   = inv cc
      aaT  = tr aa; bbT   = tr bb
      q    = ww * bbT * gg * (p - aa * x)
      ee   = - cc * aaT * gg * bb * ww
      dd   = ww + sw ee uu
      r    = p - aa*x - bb*q
      ch   = scalar $ sw r gg

      gb   = gg - sw gg (sw bbT ww)
      uu'  =  uu - sw aa gb
      duu  = det uu'
      bad  = duu < 0
      cx   = if bad then 1000.0 `debug` ("--> ksm' bad" ++ show uu')
                    else cx'' where
                      cc'  = inv uu' -- `debug` ("--> ksm " ++ show uu')
                      x'   = cc' * (uu*x - aaT * gb * p)
                      dx   = x - x'
                      cx'  = scalar $ sw dx uu'
                      cx'' = if cx' < 0 then 2000.0 `debug` ("--> ksm chi2 is " ++ show cx' ++ ", " ++ show ch ++ ", " ++ show ((max cx' 0) + ch))
                                        else cx'
      chi2 = cx + ch
  return (QMeas q dd w0, chi2)


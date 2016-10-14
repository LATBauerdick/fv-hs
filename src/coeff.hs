-- file src/coeff.hs
--
module Coeff ( w2pt, fvABh0, h2p4, q2p4, invMass, showErr ) where

import Data.Matrix ( (!), getDiag, getCol )
import Data.Vector ( zip )
import Types ( M, M33, V3, MMeas (..), HMeas (..), QMeas (..), PMeas (..)
  , ABh0 (..) )
import Matrix ( sub, sub2, toList, fromList, fromList2 )
import Text.Printf

w2pt :: Double
w2pt = 4.5451703E-03

mπ :: Double
mπ = 0.1395675E0

--fvq :: V -> V -> V
-- should calc q at nearest point of approach to v
--fvq h v = Data.Vector.take 3 h

fvABh0 :: M -> M -> ABh0
fvABh0 v h = ABh0 aa bb h0 where
  aa = fromList2 5 3 [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0,0,0, 0,0,0]
  bb = fromList2 5 3 [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0,0,0, 0,0,0]
  h0 = h

-- print a imomentum vector with error
showErr :: String -> PMeas -> IO ()
showErr s (PMeas p cp) = do
  putStr s
  let
    s2p        = getDiag cp
    f (x, s2)  = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
      where dx = sqrt s2
    in mapM_ f $ Data.Vector.zip (getCol 1 p) s2p
  putStrLn " GeV"

q2p4 :: QMeas -> PMeas
q2p4 (QMeas q cq) = h3p q cq

h2p4 :: HMeas -> PMeas
h2p4 (HMeas h ch) =
  h3p (sub 3 h) (sub2 3 ch)

h3p :: V3 -> M33 -> PMeas
h3p h3 ch = (PMeas p0 cp0) where
  m = mπ
  [w,tl,psi0] = toList 3 h3
  sph  = sin psi0
  cph  = cos psi0
  pt   = w2pt / abs w
  px   = pt * cph
  py   = pt * sph
  pz   = pt * tl
  e = sqrt(px^2 + py^2 + pz^2 + m^2)
  ps = w2pt / w
  dpdk = ps*ps/w2pt
  xy = 2.0*ps*dpdk*cph*sph*ch!(1,3)
  sxx = (dpdk*cph)^2 * ch!(1,1) + (ps*sph)^2 * ch!(3,3) + xy
  sxy = cph*sph*(dpdk*dpdk*ch!(1,1) - ps*ps*ch!(3,3)) +
           ps*dpdk*(sph*sph-cph*cph)*ch!(1,3)
  syy = (dpdk*sph)^2 * ch!(1,1) + (ps*cph)^2 * ch!(3,3) - xy
  sxz = dpdk*dpdk*cph*tl*ch!(1,1) -
           ps*dpdk*(cph*ch!(1,2)-sph*tl*ch!(1,3)) -
           ps*ps*sph*ch!(2,3)
  syz = dpdk*dpdk*sph*tl*ch!(1,1) -
           ps*dpdk*(sph*ch!(1,2) + cph*tl*ch!(1,3)) +
           ps*ps*cph*ch!(2,3)
  szz = (dpdk*tl)^2 * ch!(1,1) + ps*ps*ch!(2,2) -
           2.0*ps*dpdk*tl*ch!(1,2)
  sxe = (px*sxx + py*sxy + pz*sxz)/e
  sye = (px*sxy + py*syy + pz*syz)/e
  sze = (px*sxz + py*syz + pz*szz)/e
  see = (px*px*sxx + py*py*syy + pz*pz*szz +
         2.0*(px*(py*sxy + pz*sxz) + py*pz*syz))/e/e

  p0 = fromList 4 [px,py,pz,e]
  cp0 = fromList2 4 4 [sxx, sxy, sxz, sxe, sxy, syy, syz, sye, sxz, syz, szz, sze, sxe, sye, sze, see]

invMass :: [PMeas] -> MMeas
invMass pl@(h:t) = mass ptot where
  ptot = foldr sumP h t


sumP :: PMeas -> PMeas -> PMeas
sumP (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)

mass :: PMeas -> MMeas
mass (PMeas p cp) = mm  where
  [px,py,pz,e] = toList 4 p
  c11              = cp!(1,1)
  c12              = cp!(1,2)
  c13              = cp!(1,3)
  c14              = cp!(1,4)
  c22              = cp!(2,2)
  c23              = cp!(2,3)
  c24              = cp!(2,4)
  c33              = cp!(3,3)
  c34              = cp!(3,4)
  c44              = cp!(4,4)
  m                = sqrt $ max (e^2-px^2-py^2-pz^2) 0
  sigm0            = px*c11*px + py*c22*py + pz*c33*pz + e*c44*e +
                       2.0*(px*(c12*py + c13*pz - c14*e)
                          + py*(c23*pz - c24*e)
                          - pz*c34*e)
  sigm             = ( sqrt $ max sigm0 0 ) / m
  mm               = MMeas m sigm

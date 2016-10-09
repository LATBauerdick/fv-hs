-- file src/coeff.hs
--
module Coeff ( w2pt, h2p4, q2p4, invMass ) where

import Data.Matrix ( identity, fromLists, elementwise, (!), submatrix )
import Data.Vector ( take )
import qualified Data.Vector ( fromList, (!) )
import Types ( M33, V3, MMeas (..), HVec (..), QVec (..), PMeas (..) )

w2pt :: Double
w2pt = 4.5451703E-03

mπ :: Double
mπ = 0.1395675E0

q2p4 :: QVec -> PMeas
q2p4 (QVec q cq) = h3p q cq

h2p4 :: HVec -> PMeas
h2p4 (HVec h ch) =
  h3p (Data.Vector.take 3 h) (submatrix 1 3 1 3 ch)

h3p :: V3 -> M33 -> PMeas
h3p h3 ch = (PMeas p0 cp0) where
  m = mπ
  w    = h3 Data.Vector.! 0
  tl   = h3 Data.Vector.! 1
  psi0 = h3 Data.Vector.! 2
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

  p0 = Data.Vector.fromList [px,py,pz,e ]
  cp0 = fromLists [[sxx, sxy, sxz, sxe], [sxy, syy, syz, sye], [sxz, syz, szz, sze], [sxe, sye, sze, see]]

invMass :: [PMeas] -> MMeas
invMass pl@(h:t) = mass ptot where
  ptot = foldr sumP h t


sumP :: PMeas -> PMeas -> PMeas
sumP (PMeas p cp) (PMeas p' cp') = PMeas (Data.Vector.fromList [p0, p1, p2, p3]) sumA where
  p0   = (p Data.Vector.! 0) + (p' Data.Vector.! 0)
  p1   = (p Data.Vector.! 1) + (p' Data.Vector.! 1)
  p2   = (p Data.Vector.! 2) + (p' Data.Vector.! 2)
  p3   = (p Data.Vector.! 3) + (p' Data.Vector.! 3)
  sumA = elementwise (+) cp cp'     -- sum up the covariance matrices

mass :: PMeas -> MMeas
mass (PMeas p cp) = mm  where
  px               = p Data.Vector.! 0
  py               = p Data.Vector.! 1
  pz               = p Data.Vector.! 2
  e                = p Data.Vector.! 3
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

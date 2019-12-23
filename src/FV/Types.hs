-- {-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE ExplicitForAll #-}
--{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
-- {-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
--{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RankNTypes #-}
--{-# LANGUAGE RebindableSyntax #-}
--{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE OverloadedLists #-}
--{-# LANGUAGE NamedFieldPuns #-}

module FV.Types
  ( MCtruth (..)
  , Prong (..)
  , Chi2 (..)
  , VHMeas (..), hFilter
  , XMeas (..), vBlowup, xXMeas, yXMeas, zXMeas
  , DMeas (..), Pos, distance
  , HMeas (..)
  , QMeas (..), fromHMeas
  , PMeas (..), fromQMeas, invMass
  , MMeas (..)
  , XFit (..)
    , chi2Vertex, zVertex, dzVertex, rVertex, drVertex, z0Helix
  ) where

import Prelude.Extended
import qualified Data.Vector.Unboxed as A
  ( drop, zip, map, foldl )
import Data.Foldable ( fold )

import Data.Cov
import Data.Cov.Jac ( Jac (..) )

-----------------------------------------------
-- Prong
-- a prong results from a vertex fit of N helices
newtype Chi2  = Chi2 Number
instance Show Chi2 where
  show (Chi2 c2) = show c2
instance Semiring Chi2 where
  add (Chi2 c0) (Chi2 c1) = Chi2 (c0+c1)
  zero = Chi2 0.0
  mul (Chi2 c0) (Chi2 c1) = Chi2 (c0*c1)
  one = Chi2 1.0
instance Num Chi2 where
  fromInteger i = Chi2 $ fromInteger i
  negate (Chi2 c) = Chi2 (negate c)
  abs (Chi2 c) = Chi2 (abs c)
  signum (Chi2 c) = Chi2 (signum c)
  (*) = error "cannot multiply Chi2*Chi2 to return a Chi2"
  (+) (Chi2 a) (Chi2 b) = Chi2 (a+b)
data Prong = Prong
          { nProng        :: Int
          , fitVertex     :: XMeas
          , fitMomenta    :: List QMeas
          , fitChi2s      :: List Chi2
          , measurements  :: VHMeas
          }
-- fitMomenta :: Prong -> List QMeas
-- fitMomenta (Prong {fitMomenta= f}) = f
instance Show Prong where
  show _ = "Prong!!!!!!!!"

-----------------------------------------------
-- VHMeas
--
data VHMeas = VHMeas {
    vertex      :: XMeas
  , helices     :: List HMeas
                     }
-- vertex :: VHMeas -> XMeas
-- vertex (VHMeas {vertex=v}) = v
-- helices :: VHMeas -> List HMeas
-- helices (VHMeas {helices=hs}) = hs
{-- instance Semigroup VHMeas where --}
{--   append (VHMeas v hs) (VHMeas _ hs') = VHMeas v ( hs ++ hs' ) --}
{--   mempty = VHMeas (XMeas (Matrix.zero 3 1) (Matrix.zero 3 3)) [] --}

instance Show VHMeas where
  show VHMeas {vertex=v_, helices=hs} = "VHMeas w/ " <> show (length hs)
                                    <> " tracks. " <> show v_

vBlowup :: Number -> VHMeas -> VHMeas
{-- vBlowup scale vm = over vertexLens (blowup scale) vm where --}
vBlowup scale VHMeas {vertex= v_, helices= hs} =
    VHMeas {vertex= blowup scale v_, helices= hs} where
  blowup :: Number -> XMeas -> XMeas -- blow up diag of cov matrix
  blowup s (XMeas v__ cv) = XMeas v__ cv' where
    cv' = scaleDiag s cv

hFilter :: List Int -> VHMeas -> VHMeas
hFilter is VHMeas {vertex=v_, helices=hs} = VHMeas {vertex=v_, helices= iflt is hs}

hRemove :: Int -> VHMeas -> VHMeas
hRemove indx VHMeas {vertex=v_, helices=hs} = VHMeas {vertex=v_, helices=irem indx hs}

-----------------------------------------------
-- MCtruth
--
newtype MCtruth = MCtruth {
    pu_zpositions :: List Number
                       }
instance Show MCtruth where
  show MCtruth {pu_zpositions=puz} = "MCtruth w/" <> show (length puz)
                                                <> " PU z positions."


-----------------------------------------------
-- HMeas
-- 5-vector and covariance matrix for helix measurement
--
data HMeas = HMeas Vec5 Cov5 Number
instance Show HMeas where
  show (HMeas h ch _) = s' where
    sh = A.map sqrt $ diag ch
    hs = toArray h
    s00 = to5fix x <> " +-" <> to5fix dx where
      x  = uidx hs 0
      dx = uidx sh 0
    s' = A.foldl f s00 (A.drop 1 $ A.zip hs sh) where
      f s (x, dx)  = s <> to3fix x <> " +-" <> to3fix dx

z0Helix :: HMeas -> Number
z0Helix (HMeas h _ _) = z0 where
  v = val h
  z0 = uidx v 4
-----------------------------------------------
-- QMeas
-- 3-vector and covariance matrix for momentum measurement
--
mπ :: Number
mπ = 0.1395675
data QMeas = QMeas Vec3 Cov3 Number
instance Show QMeas where
  show = showQMeas
-- print QMeas as a 4-momentum vector with errors, use pt and pz
showQMeas :: QMeas -> String
showQMeas (QMeas q cq w2pt) = s' where
  f :: String -> (Number, Number) -> String
  f s (x, dx)  = s <> to3fix x <> " +-" <> to3fix dx
  m          = mπ
  wp         = w2pt
  qs :: Array Number
  qs         = toArray q
  w          = uidx qs 0
  tl         = uidx qs 1
  psi0       = uidx qs 2
  pt         = wp / abs w
  pz         = pt*tl
  psi        = psi0*180.0/pi
  e          = sqrt(pt*pt  + pz*pz + m*m)
  jj :: Jac34
  jj         = Jac { v = [ -wp/w/w, -wp/w/w*tl, 0.0, -(pz*pz + pt*pt)/w/e
                          , 0.0, wp/w, 0.0, pt*pt*tl/e
                          , 0.0, 0.0, 1.0, 0.0], nr = 3 }
  cq'        = jj .*. cq
  p'         = [pt, pz, psi, e]
  dp         = A.map sqrt $ diag cq'
  d1         = uidx dp 0
  d2         = uidx dp 1
  d3         = uidx dp 2
  d4         = uidx dp 3
  dp'        = [d1, d2, d3*180.0/pi, d4]
  s'         = A.foldl f "" ( A.zip p' dp' ) <> " GeV"

fromHMeas :: HMeas -> QMeas -- just drop the d0, z0 part... fix!!!!
fromHMeas (HMeas h ch w2pt) = QMeas q cq w2pt where
  q = subm 3 h
  cq = subm2 3 ch

h2p :: HMeas -> PMeas
h2p = fromQMeas <<< fromHMeas

-----------------------------------------------
-- PMeas
-- 4-vector and coavariance matrix for momentum px,py,pz and energy
--
data PMeas = PMeas Vec4 Cov4
instance Semigroup PMeas where
  (<>) (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)
instance Monoid PMeas where
  mappend (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)
  mempty = PMeas (fromArray [0.0,0.0,0.0,0.0]) zero
instance Show PMeas where
  show = showPMeas
-- print PMeas as a 4-momentum vector px,py,pz,E with errors
showPMeas :: PMeas -> String
showPMeas (PMeas p cp) = s' where
  sp         = A.map sqrt $ diag cp
  f s (x, dx)  = s <> to3fix x <> " +-" <> to3fix dx -- \xc2b1 ±±±±±
  s' = A.foldl f "" (A.zip (toArray p) sp) <> " GeV"

invMass :: List PMeas -> MMeas
invMass = pmass <<< fold

pmass :: PMeas -> MMeas
pmass (PMeas p cp) = mm  where
  ps    = toArray p
  px    = uidx ps 0
  py    = uidx ps 1
  pz    = uidx ps 2
  e     = uidx ps 3
  cps   = toArray cp
  c11   = uidx cps 0
  c12   = uidx cps 1
  c13   = uidx cps 2
  c14   = uidx cps 3
  c22   = uidx cps 5
  c23   = uidx cps 6
  c24   = uidx cps 7
  c33   = uidx cps 10
  c34   = uidx cps 11
  c44   = uidx cps 15
  m     = sqrt $ max (e*e-px*px-py*py-pz*pz) 0.0
  sigm0 = px*c11*px + py*c22*py + pz*c33*pz + e*c44*e +
            2.0*(px*(c12*py + c13*pz - c14*e)
               + py*(c23*pz - c24*e)
               - pz*c34*e)
  dm    =  sqrt ( max sigm0 0.0 ) / m
  mm    = MMeas {m=m, dm=dm}

fromQMeas :: QMeas -> PMeas
fromQMeas (QMeas q0 cq0 w2pt) = PMeas p0 cp0 where
  m = mπ
  q0s = toArray q0
  w    = uidx q0s 0
  tl   = uidx q0s 1
  psi0 = uidx q0s 2
  sph  = sin psi0
  cph  = cos psi0
  pt   = w2pt / abs w
  px   = pt * cph
  py   = pt * sph
  pz   = pt * tl
  sqr x= x*x
  e    = sqrt(px*px + py*py + pz*pz + m*m)
  ps   = w2pt / w
  dpdk = ps*ps/w2pt
  cq0s = toArray cq0
  c11  = uidx cq0s 0
  c12  = uidx cq0s 1
  c13  = uidx cq0s 2
  c22  = uidx cq0s 4
  c23  = uidx cq0s 5
  c33  = uidx cq0s 8
  xy   = 2.0*ps*dpdk*cph*sph*c13
  sxx  = sqr (dpdk*cph) * c11 + sqr (ps*sph) * c33 + xy
  sxy  = cph*sph*(dpdk*dpdk*c11 - ps*ps*c33) +
           ps*dpdk*(sph*sph-cph*cph)*c13
  syy  = sqr (dpdk*sph) * c11 + sqr (ps*cph) * c33 - xy
  sxz  = dpdk*dpdk*cph*tl*c11 -
           ps*dpdk*(cph*c12-sph*tl*c13) -
           ps*ps*sph*c23
  syz  = dpdk*dpdk*sph*tl*c11 -
           ps*dpdk*(sph*c12 + cph*tl*c13) +
           ps*ps*cph*c23
  szz  = sqr (dpdk*tl) * c11 + ps*ps*c22 -
           2.0*ps*dpdk*tl*c12
  sxe  = (px*sxx + py*sxy + pz*sxz)/e
  sye  = (px*sxy + py*syy + pz*syz)/e
  sze  = (px*sxz + py*syz + pz*szz)/e
  see  = (px*px*sxx + py*py*syy + pz*pz*szz +
         2.0*(px*(py*sxy + pz*sxz) + py*pz*syz))/e/e

  cp0  = fromArray [ sxx, sxy, sxz, sxe
                        , syy, syz, sye
                             , szz, sze
                                  , see]
  p0   = fromArray [px,py,pz,e]

-----------------------------------------------
-- MMeas
-- scalar mass and error calculated from PMeas covariance matrices
--
data MMeas = MMeas
            { m :: Number
            , dm :: Number
            }
instance Show MMeas where
  show MMeas {m=m, dm=dm} = " " <> to1fix (m*1000.0) <> " +-" <> to1fix (dm*1000.0) <> " MeV"

-----------------------------------------------
-- XMeas
-- 3-vector and covariance matrix for position/vertex measurement
--
data DMeas = DMeas Number Number -- distance and error
instance Show DMeas where
  show (DMeas d sd) = to2fix d <> " +-" <> to2fix sd <> " cm"
class Pos a where
  distance :: a -> a -> DMeas -- distance between two positions
instance Pos XMeas where
  distance (XMeas v0 vv0) (XMeas v1 vv1) = DMeas d sd where
    v0s = toArray v0
    x0 = uidx v0s 0
    y0 = uidx v0s 1
    z0 = uidx v0s 2
    v1s = toArray v1
    x1 = uidx v1s 0
    y1 = uidx v1s 1
    z1 = uidx v1s 2

    d  = sqrt(sqr(x0-x1) + sqr(y0-y1) + sqr(z0-z1))
    dd :: Vec3
    dd = fromArray [(x0-x1)/d, (y0-y1)/d, (z0-z1)/d]
    tem0 = dd .*. vv0
    tem1 = dd .*. vv1
    sd   = sqrt (tem0 + tem1)

data XMeas = XMeas Vec3 Cov3
instance Semigroup XMeas where
  (<>) (XMeas x1 cx1) (XMeas x2 cx2) = XMeas (x1 + x2) (cx1 + cx2)
instance Monoid XMeas where
  mempty = XMeas zero zero
  mappend = undefined
instance Show XMeas where
  show = showXMeas
-- return a string showing vertex position vector with errors
showXMeas :: XMeas -> String
showXMeas (XMeas v cv) = s' where
  vv         = toArray v
  x          = uidx vv 0
  y          = uidx vv 1
  z          = uidx vv 2
  s2v        = A.map sqrt $ diag cv
  dx         = uidx s2v 0
  dy         = uidx s2v 1
  dz         = uidx s2v 2
  f :: Number -> Number -> String -> String
  f x dx s  = s <> to2fix x <>  " +-" <> to2fix dx
  s' = f z dz <<< f y dy <<< f x dx $
    "(r,z) =" <> "(" <> to2fix (sqrt (x*x + y*y))
              <> ", " <> to2fix z <> "), x y z ="

xXMeas :: XMeas -> Number
xXMeas (XMeas v _) = x where
  x = uidx (val v) 0
yXMeas :: XMeas -> Number
yXMeas (XMeas v _) = y where
  y = uidx (val v) 1
zXMeas :: XMeas -> Number
zXMeas (XMeas v _) = z where
  z = uidx (val v) 2

data XFit = XFit Vec3 Cov3 Chi2
instance Show XFit where
  show (XFit x xx c2) = showXMeas (XMeas x xx) <> ", chi2=" <> show c2

chi2Vertex :: XFit -> Chi2
chi2Vertex (XFit _ _ c2) = c2

zVertex :: XFit -> Double
zVertex (XFit v _ _) = z where
  z = uidx (val v) 2

dzVertex :: XFit -> Double
dzVertex (XFit _ dv _) = dz where
  dz = sqrt $ uidx (val dv) 5

rVertex :: XFit -> Double
rVertex (XFit v _ _) = r where
  r = sqrt $ (sqr (uidx (val v) 0)) + (sqr (uidx (val v) 1))

drVertex :: XFit -> Double
drVertex (XFit v _ _) = const 0.002 dr where
  dr = sqrt $ (sqr (uidx (val v) 0)) + (sqr (uidx (val v) 1))

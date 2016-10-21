-- file src/Fit.hs
module Types (
  M, V, M33, V3, V5, C44, XMeas (..), HMeas (..), QMeas (..)
             , PMeas (..), MMeas (..), Prong (..), VHMeas (..)
             , ABh0 (..), Chi2
             , showXMeas, showPMeas, showQMeas, showMMeas
             , w2pt, mπ
             ) where

import Text.Printf
import qualified Data.Matrix ( Matrix, getDiag, getCol, toList, zero
                             , fromLists, transpose )
import qualified Data.Vector ( zip, map, fromList, toList )

w2pt :: Double
w2pt = 4.5451703E-03

mπ :: Double
mπ = 0.1395675E0


type M     = Data.Matrix.Matrix Double
type V     = Data.Matrix.Matrix Double
type M33   = Data.Matrix.Matrix Double
type V3    = Data.Matrix.Matrix Double
type V5    = Data.Matrix.Matrix Double
type N     = Int
type Chi2  = Double
data Prong = Prong N XMeas [QMeas] [Chi2] deriving Show -- a prong results from a vertex fit of N helices

data VHMeas a = VHMeas XMeas [a] deriving Show
instance Monoid (VHMeas a) where
  mappend (VHMeas v1 as1) (VHMeas v2 as2) = VHMeas (v1) ( as1 ++ as2 )
instance Functor VHMeas where -- apply f to each a
  fmap f (VHMeas v (a:as)) = VHMeas v ((f a):(fmap f as))

data ABh0 = ABh0 M M M

type X3 = V
type C33 = M
data XMeas = XMeas X3 C33 deriving Show -- 3-vector and covariance matrix for position/vertex measurement

type H5 = V
type C55 = M
data HMeas = HMeas H5 C55 deriving Show -- 5-vector and covariance matrix for helix measurement

type Q3 = V
data QMeas = QMeas Q3 C33 deriving Show -- 3-vector and covariance matrix for momentum measurement

-- 4-vector and coavariance matrix for momentum px,py,pz and energy
type P4 = V -- four-vector
type C44 = M -- 4x4 covariance matrix
data PMeas = PMeas P4 C44 deriving Show
instance Monoid (PMeas) where
  mappend (PMeas p1 cp1) (PMeas p2 cp2) = PMeas (p1+p2) (cp1 + cp2)
  mempty = PMeas (Data.Matrix.fromLists [[0.0,0.0,0.0,0.0]])
            ((Data.Matrix.zero 4 4)::C44)
-- instance Functor PMeas where
--   fmap f (PMeas p cp) = f p cp

type D = Double
data MMeas = MMeas D D deriving Show -- mass and error

-- show instances -- refactor!!

-- print a vertext position vector with errors
showXMeas :: String -> XMeas -> IO ()
showXMeas s (XMeas v cv) = do
  putStr s
  let
    s2v        = Data.Vector.map sqrt $ Data.Matrix.getDiag cv
    f (x, dx)  = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
    in mapM_ f $ Data.Vector.zip (Data.Matrix.getCol 1 v) s2v
  putStrLn " cm"

showMMeas :: String -> MMeas -> IO ()
showMMeas s (MMeas m dm) = do
  putStr s
  printf "%8.3f ± %8.3f" (m::Double) (dm::Double)
  putStrLn " GeV"

-- print PMeas as a 4-momentum vector px,py,pz,E with errors
showPMeas :: String -> PMeas -> IO ()
showPMeas s (PMeas p cp) = do
  putStr s
  let
    sp         = Data.Vector.map sqrt $ Data.Matrix.getDiag cp
    f (x, dx)  = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
    in mapM_ f $ Data.Vector.zip (Data.Matrix.getCol 1 p) sp
  putStrLn " GeV"


-- print QMeas as a 4-momentum vector with errors, use pt and pz
showQMeas :: String -> QMeas -> IO ()
showQMeas s (QMeas q cq) = do
  putStr s
  let
    f (x, dx)    = printf "%8.3f ± %8.3f" (x::Double) (dx::Double)
    m = mπ
    wp = w2pt
    [w,tl,psi0] = take 3 (Data.Matrix.toList q)
    pt   = wp / abs w
    pz = pt*tl
    psi = psi0*180.0/pi
    e = sqrt(pt^2  + pz^2 + m^2)
    jj   = Data.Matrix.fromLists [
            [-wp/w/w, -wp/w/w*tl,0, -(pz*pz + pt*pt)/w/e ]
          , [0, wp/w, 0, pt*pt*tl/e]
          , [0, 0, 1.0, 0] ]
    cq'  = (Data.Matrix.transpose jj) * cq* jj
    p'   = Data.Vector.fromList [pt, pz, psi, e]
--    dp'  = Data.Vector.map sqrt $ Data.Matrix.getDiag cq'
    [d1,d2,d3,d4]  = Data.Vector.toList $ Data.Vector.map sqrt $ Data.Matrix.getDiag cq'
    dp' = Data.Vector.fromList [d1,d2,d3*180.0/pi,d4]
    in mapM_ f $ Data.Vector.zip p' dp'
  putStrLn " GeV"


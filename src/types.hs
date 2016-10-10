-- file src/Fit.hs
module Types (
              M, V, M33, V3, XVec (..), HVec (..), QVec (..)
            , PMeas (..), MMeas (..), Prong (..)
             ) where

import Data.Matrix ( Matrix )
import Data.Vector ( Vector )

type M = Matrix Double
type V = Vector Double
type M33 = Matrix Double
type V3 = Vector Double
type N = Int
type Chi2 = Double
data Prong = Prong N XVec [QVec] [Chi2] deriving Show -- a prong results from a vertex fit of N helices


type X3 = V
type C33 = M
data XVec = XVec X3 C33 deriving Show -- 3-vector and covariance matrix for position/vertex measurement

type H5 = V
type C55 = M
data HVec = HVec H5 C55 deriving Show -- 5-vector and covariance matrix for helix measurement

type Q3 = V
data QVec = QVec Q3 C33 deriving Show -- 3-vector and covariance matrix for momentum measurement

type P4 = V -- four-vector
type C44 = M -- 4x4 covariance matrix
data PMeas = PMeas P4 C44 deriving Show -- 4-vector and coavariance matrix for momentum px,py,pz and energy

type D = Double
data MMeas = MMeas D D deriving Show -- mass and error

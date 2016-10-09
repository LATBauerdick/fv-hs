-- file src/Fit.hs
module Fit ( fit ) where

import Data.Matrix ( Matrix, submatrix )
import Data.Vector ( take )
import Types ( XVec (..), HVec (..), QVec (..), Prong (..) )

fvq :: HVec -> QVec
fvq (HVec h ch) = q where
  q = QVec (Data.Vector.take 3 h) (submatrix 1 3 1 3 ch)

fit :: XVec -> [HVec] -> Prong
fit v0 hl = pr where
  pr = Prong 6 v0 (map fvq hl) [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]


module Test.Probability ( doProbability ) where

import Prelude
import Control.Parallel.Strategies
import System.Random
import Data.Random.Normal ( normals )
import Statistics.Sample ( meanVariance )

import qualified Data.Vector.Unboxed as V ( Vector, fromList, toList )

import Data.List ( foldl', unfoldr, mapAccumL, (!!) )

import qualified Numeric.Probability.Distribution as Dist
import qualified Numeric.Probability.Transition as Trans
import qualified Numeric.Probability.Trace as Trace
import qualified Numeric.Probability.Random as Rnd
import Numeric.Probability.Simulation ((~..), (~*.), )
import Numeric.Probability.Percentage
    (Dist, Trans, RTrans, Expand, RExpand, Space, )

import FV.Types ( XMeas (..), VHMeas (..), HMeas (..), Prong (..), MMeas (..)
             , invMass, q2p, v3,l3,v5,l5 )
import Data.SimpleMatrix ( toList, fromList, chol, scalar )
import FV.Fit ( fit, fitw )

import qualified Graphics.Gnuplot.Frame.OptionSet as Opts
import Graphics.Histogram

type Height = Int

data Tree = Alive Height | Hit Height | Fallen
            deriving (Ord,Eq,Show)

grow :: Trans Tree
grow (Alive h) = Dist.normal (map Alive [h+1..h+5])
grow _ = error "TreeGrowth.grow: only alive trees can grow"

hit :: Trans Tree
hit (Alive h) = Dist.certainly (Hit h)
hit _ = error "TreeGrowth.hit: only alive trees can be hit"

fall :: Trans Tree
fall _ = Dist.certainly Fallen

evolve :: Trans Tree
evolve t =
   case t of
      (Alive _) -> Trans.unfold (Dist.enum [90,4,6] [grow,hit,fall]) t
--    (Alive _) -> Trans.unfold (Dist.relative [0.9,0.04,0.06] [grow,hit,fall]) t
      _         -> Dist.certainly t

{- |
tree growth simulation:
 start with seed and run for n generations
-}
seed :: Tree
seed = Alive 0

-- * exact results

-- | @hist n@ : history of tree distributions for n generations
hist :: Int -> Expand Tree
hist n = Trace.walk n (evolve =<<) . return

simHist :: Int -> Int -> RExpand Tree
simHist k n = (k,n) ~.. evolve
-- simTree k n = k ~. n *. random evolve
-- simTree k n = (k,n) ~*. evolve

h2 :: Space Tree
h2  = hist 2 seed

sh2 :: IO ()
sh2 = Rnd.print $ simHist 2000 5 seed


doProbability :: Int -> VHMeas -> IO ()
doProbability cnt vm = do
  putStrLn $ "----------------------------------------------"
  print $ hist 2 seed
  putStrLn $ "----------------------------------------------"
  Rnd.print $ simHist 2000 2 seed
  putStrLn $ "----------------------------------------------"


  putStrLn $ "Fit Mass  " ++ (show . invMass . map q2p . fitMomenta . fit $ vm)


  g <- newStdGen
  let hf :: V.Vector Double
      hf = V.fromList
          . withStrategy (parBuffer 100 rdeepseq)
          . map fitm
          . take cnt . gen vm . normals $ g -- produce a list of cnt randomized VHMeas
      (mean, var) = meanVariance hf
  putStrLn $ "Mean Mass " ++ show (MMeas mean (sqrt var))
  let hist = histogram binSturges (V.toList hf)
  _ <- plot "invMass.png" hist
  putStrLn $ "----------------------------------------------"
  return ()

gen :: VHMeas -> [Double] -> [VHMeas]
gen v rs = v' : gen v rs' where
  (v', rs') = randVH v rs

-- calc fitted invariant mass of VHMeas
fitm :: VHMeas -> Double
fitm vm = m where
  (MMeas m _) = invMass . map q2p . fitMomenta $ fit vm

-- randomize the helices in the supplied VHMeas
-- and return randomized VHMeas and remaining randoms list
randVH :: VHMeas -> [Double] -> (VHMeas, [Double])
randVH (VHMeas v hl) rs = (VHMeas v hl', rs') where
  (rs', hl') = mapAccumL randH rs hl

-- randomize a single helix parameters measurement, based on the cov matrix
-- return randomized helix and "remaining" random numbers
randH :: [Double] -> HMeas -> ([Double], HMeas)
randH (r0:r1:r2:r3:r4:rs) (HMeas h hh w0) = (rs, HMeas h' hh w0) where
  h' = v5 $ zipWith (+) (l5 h) (l5 (chol hh * v5 [r0,r1,r2,r3,r4]))


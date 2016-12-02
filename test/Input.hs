-- file: test/Input.hs
module Input ( hSlurp, hSlurpAll, dataFiles ) where

import Types ( HMeas (..), XMeas (..), VHMeas (..) )
import Matrix ( fromList, fromList2, scaleDiag, sw )
import System.Directory
import qualified Data.List.Split ( chunksOf )
import Data.Monoid
import Debug.Trace ( trace )
debug :: a -> String -> a
debug = flip trace

-- slurp in the measurements of vertex and helices
-- from a list of Doubles
hSlurp' :: [Double] -> VHMeas HMeas
hSlurp' inp = VHMeas v hl where
  v0        = fromList 3 $ take 3 inp   -- initial vertex pos
  cv0       = scaleDiag 10000.0 $ fromList2 3 3 $ take 9 $ drop 3 inp -- cov matrix
  v         = XMeas v0 cv0
  w2pt      = inp !! 12                 -- how to calc pt from w; 1 in case of CMS
  nt        = round (inp !! 13) ::Int    -- number of helices to follow
  hl        = map nxtH . Data.List.Split.chunksOf 30 $ drop 14 inp -- list of helix params and cov
  nxtH :: [Double] -> HMeas
  nxtH ds =  hm where
    (ih, ich) = splitAt 5 ds
    hm = if w2pt /= 1.0
            then  let
                    h'   = fromList 5 ih
                    ch'  = fromList2 5 5 ich
                    in HMeas h' ch' w2pt
            else let
              -- FV works in terms of a perigee system
              -- w = omega = 1/R is curvature radius
              -- tl = tan lambda = tangent of dipping angle (0 for pt-max)
              -- psi = angle in xy plane
              -- d0 = distance of closest approach to 0 in xy plan
              -- z0 = positiion on z axis
              --
              -- CMS uses instead
              -- q/p = charge over momentum
              -- theta = dipping angle
              -- etc
              --
              w0               = 0.003*3.8  -- CMS case: field is 3.8 T, give R in cm
              [h0,h1,h2,h3,h4] = ih
              st               = sin h1
              ct               = cos h1
              w                = h0 * w0 / ct
              tl               = st / ct
              pt               = w0/w
              j00              = w0 / ct -- `debug` ((show (pt*tl)) )
              j01              = h0 * w0 * st/ct/ct
              j11              = 1.0 / ct / ct
              j10              = 0
              jj               = fromList2 5 5 [  j00, j01, 0, 0, 0
                                                , j10, j11, 0, 0, 0
                                                , 0, 0, 1.0, 0, 0
                                                , 0, 0, 0, 1.0, 0
                                                , 0, 0, 0, 0, 1.0 ]
              h'               = fromList 5 (w : tl : h2 : h3 : h4 : [])
              ch'              = sw jj (fromList2 5 5 ich)
              in HMeas h' ch' w0
-- slurps up a bunch of Doubles from a text data file into a list
-- and parses them w/ hSlurp' to a vertex and a set of helix measurements
hSlurp :: String -> IO (VHMeas HMeas)
hSlurp path = do
  ds <- readFile path
  return $ hSlurp' $ map readDouble (words ds)
    where readDouble = read :: String -> Double

-- slurp all files named in a list of pathNames
hSlurpAll :: [String] -> IO (VHMeas HMeas)
hSlurpAll (x:[]) = hSlurp x
hSlurpAll (x:xs) = do
  v0 <- hSlurp x
  vs <- hSlurpAll xs
  return $ v0 <> vs

justFiles :: FilePath -> Bool
justFiles f = f /= "." && f /= ".."

dataFiles :: FilePath -> IO [FilePath]
dataFiles path = do
  fileNames <- filter justFiles <$> getDirectoryContents path
  let addHead = (++) (path ++ "/")
  return $ map addHead fileNames

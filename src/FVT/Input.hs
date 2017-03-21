-- file: test/Input.hs
module FVT.Input ( hSlurp, hSlurpAll, dataFiles ) where

import Prelude
import FV.Types ( HMeas (..), XMeas (..), VHMeas (..), MCtruth (..) )
import FV.Matrix ( fromList, fromList2, sw )
import System.Directory
import qualified Data.List.Split ( chunksOf )
import Data.Monoid
-- import Debug.Trace ( trace )
-- debug :: a -> String -> a
-- debug = flip trace

-- slurp in the measurements of vertex and helices
-- from a list of Doubles
hSlurp' :: [Double] -> VHMeas
hSlurp' inp = VHMeas v hl where
  v0        = fromList 3 $ take 3 inp   -- initial vertex pos
  cv0       = fromList2 3 3 $ take 9 $ drop 3 inp -- cov matrix
  v         = XMeas v0 cv0
  w2pt      = inp !! 12                 -- how to calc pt from w; 1 in case of CMS
  nt        = round (inp !! 13) ::Int    -- number of helices to follow
  hl        = map nxtH . Data.List.Split.chunksOf 30 . drop 14 . take (nt*30+14) $ inp -- list of helix params and cov
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
              j00              = w0 / ct -- `debug` ((show (w0/w*tl)) )
              j01              = h0 * w0 * st/ct/ct
              j11              = 1.0 / ct / ct
              j10              = 0
              jj               = fromList2 5 5 [  j00, j01, 0, 0, 0
                                                , j10, j11, 0, 0, 0
                                                , 0, 0, 1.0, 0, 0
                                                , 0, 0, 0, 1.0, 0
                                                , 0, 0, 0, 0, 1.0 ]
              h'               = fromList 5 [w, tl, h2, h3, h4]
              ch'              = sw jj (fromList2 5 5 ich)
              in HMeas h' ch' w0
-- slurps up a bunch of Doubles from a text data file into a list
-- and parses them w/ hSlurp' to a vertex and a set of helix measurements
hSlurp :: String -> IO (VHMeas, MCtruth)
hSlurp path = do
  ds <- readFile path
  let ws = words ds
      readDouble = read :: String -> Double
      readInt = read :: String -> Int
  let (mc, ws') = if head ws == "PU_zpositions:"
                      then  let npu = readInt $ ws !! 1
                                zws = take npu $ drop 2 ws
                            in (MCtruth (map readDouble zws), drop (npu+2) ws)
                      else (MCtruth [], ws)
  return (hSlurp' $ map readDouble ws', mc)

-- slurp all files named in a list of pathNames
hSlurpAll :: [String] -> IO VHMeas
hSlurpAll [x] = do
  (v0, _) <- hSlurp x
  return v0
hSlurpAll (x:xs) = do
  (v0, _) <- hSlurp x
  vs <- hSlurpAll xs
  return $ v0 <> vs

justFiles :: FilePath -> Bool
justFiles f = f /= "." && f /= ".."

dataFiles :: FilePath -> IO [FilePath]
dataFiles path = do
  fileNames <- filter justFiles <$> getDirectoryContents path
  let addHead = (++) (path ++ "/")
  return $ map addHead fileNames

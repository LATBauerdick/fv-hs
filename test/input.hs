-- file: test/input.hs
module Input ( hSlurp, hSlurpAll, dataFiles )
  where

import Types ( HMeas (..), XMeas (..), VHMeas (..) )
import Matrix ( fromList, fromList2, scaleDiag )
import Control.Applicative
import System.Directory
import Data.Monoid

-- slurp in the measurements of vertex and helices
-- from a list of Doubles
hSlurp' :: [Double] -> VHMeas HMeas
hSlurp' inp = (VHMeas v hl) where
  v0        = fromList 3 $ take 3 inp   -- initial vertex pos
  -- cv00      = fromList2 3 3 $ take 9 $ drop 3 inp -- cov matrix
  -- cv0       = scale 10000.0 cv00
  cv0       = scaleDiag 10000.0 $ fromList2 3 3 $ take 9 $ drop 3 inp -- cov matrix
  v         = XMeas v0 cv0
  w2pt      = inp !! 12                 -- how to calc pt from w
  nt        = round (inp !! 13) ::Int    -- number of helices to follow
  i0        = drop 14 inp
  lst       = (nt-1)*30
  is        = [ drop n i0 | n <-[0,30..lst]]
  hl        = [ nxtH i | i <- is ]      -- list of helix params and cov
  nxtH :: [Double] -> HMeas
  nxtH i = HMeas h ch where
      (ih, ich) = splitAt 5 i
      h         = fromList 5 ih
      ch        = fromList2 5 5 $ take 25 ich

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

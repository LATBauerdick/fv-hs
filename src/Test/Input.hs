{-# LANGUAGE EmptyDataDecls #-}
{-# LANGUAGE ExplicitForAll #-}
--{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses #-}
--{-# LANGUAGE PartialTypeSignatures #-}
{-# LANGUAGE RankNTypes #-}
--{-# LANGUAGE RebindableSyntax #-}
--{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE OverloadedLists #-}
--{-# LANGUAGE NamedFieldPuns #-}
module Test.Input ( hSlurp, hSlurpMCtruth ) where

import Prelude
import Control.Monad (guard)
import qualified Data.Vector.Unboxed as A
  ( fromList )
import Data.List as L ( head, drop, take )
--import Data.Tuple (Tuple (..), fst)
import Data.Maybe ( mapMaybe )
--import Control.Plus (empty)

import Data.Cov
import FV.Types ( MCtruth (..), VHMeas (..), XMeas (..), HMeas (..) )
import Stuff


-- type FilePath = String -- from Node.Path

{-- listToArray :: forall a. List a -> Array a --}
{-- listToArray = ?whatGoesHere --}

-- slurp in a String of measurements of PU z-positions into MCtruth
hSlurpMCtruth :: String -> Maybe MCtruth
hSlurpMCtruth ds = mc where
  ws = words ds
  npu :: Maybe Int
  npu = do
              key <- pure <<< L.head $ ws
              guard $ key == "PU_zpositions:"
              snpu <- pure <<< L.head <<< L.drop 1 $ ws
              intFromString snpu
  mc = case npu of
              Nothing -> Nothing
              Just _ -> let
                            mcArr = mapMaybe numberFromString <<< L.drop 2 $ ws
                        in Just $ MCtruth { pu_zpositions= mcArr }

-- slurps up a String with a bunch of Doubles
-- and parses them w/ to a vertex and a set of helix measurements
hSlurp :: String -> Maybe VHMeas
hSlurp ds = vhm where
  ws = words ds
  npu :: Maybe Int
  npu = do
              key <- pure <<< L.head $ ws
              guard $ key == "PU_zpositions:"
              snpu <- pure <<< L.head <<< L.drop 1 $ ws
              intFromString snpu
  vhm = case npu of
          Nothing -> hSlurp' $ mapMaybe numberFromString ws
          Just n  -> let vArr = mapMaybe numberFromString (L.drop (n+2) ws)
                        in hSlurp' vArr

-- slurp in the measurements of vertex and helices
hSlurp' :: List Number -> Maybe VHMeas
hSlurp' inp = do
  let v0    = fromArray <<< A.fromList <<< L.take 3 $ inp       -- initial vertex pos
      cv0   = fromArray <<< A.fromList <<< L.take 9 <<< L.drop 3 $ inp -- cov matrix
      v     = XMeas v0 cv0
  w2pt      <- pure <<< L.head <<< L.drop 12 $ inp  -- how to calc pt from w; 1 in case of CMS
  mnt       <- pure <<< L.head <<< L.drop 13 $ inp  -- number of helices to follow --}
  let nt    = round mnt
      f     = case w2pt of
                  1.0 -> nxtH'        -- CMS case
                  _   -> nxtH   -- Aleph case
      hl    = mapMaybe (\i -> f w2pt (L.take 30 <<< L.drop (i*30+14) $ inp)) $ [0 .. (nt-1)]

  pure $ VHMeas { vertex= v, helices= hl }

-- get the next helix, Aleph case
nxtH :: Number -> List Number -> Maybe HMeas
nxtH w0 ds = do
  h    <- pure <<< fromArray <<< A.fromList <<< L.take 5 $ ds
  ch   <- pure <<< fromArray <<< A.fromList <<< L.take 25 <<< L.drop 5 $ ds
  pure $ HMeas h ch w0

-- get the next helix, CMS case
nxtH' :: Number -> List Number -> Maybe HMeas
nxtH' _ ds = do
  -- FV works in terms of a perigee system
  -- w = omega = 1/R is curvature radius
  -- tl = tan lambda = tangent of dipping angle (0 for pt-max)
  -- psi = angle in xy plane
  -- d0 = distance of closest approach to 0 in xy plan
  -- z0 = positiion on z axis
  --
  -- CMS uses instead
  -- q/p = charge over momentum
  -- theta = dip angle
  -- etc
  --
  [h0,h1,h2,h3,h4] <- pure <<< A.fromList <<< L.take 5 $ ds
  ch               <- pure <<< A.fromList <<< L.take 25 <<< L.drop 5 $ ds
  let w0                = 0.003*3.8  -- CMS case: field is 3.8 T, give R in cm
      st                = sin h1
      ct                = cos h1
      w                 = h0 * w0 / ct
      tl                = st / ct
      j00               = w0 / ct
      j01               = h0 * w0 * st/ct/ct
      j11               = 1.0 / ct / ct
      j10               = 0.0
      jj :: Jac55
      jj                = Jac { vj= [ j00, j01, 0.0, 0.0, 0.0
                                    , j10, j11, 0.0, 0.0, 0.0
                                    , 0.0, 0.0, 1.0, 0.0, 0.0
                                    , 0.0, 0.0, 0.0, 1.0, 0.0
                                    , 0.0, 0.0, 0.0, 0.0, 1.0
                                    ], nr= 5 }
      h' :: Vec5
      h'                = fromArray [w, tl, h2, h3, h4]
      ch' :: Cov5
      ch'               = fromArray ch
      ch''              = jj .*. ch'

  pure $ HMeas h' ch'' w0

-- slurp all files named in a list of pathNames
{-- hSlurpAll :: forall eff. --}
{--              Array FilePath --}
{--              -> Eff (fs :: FS, exception :: EXCEPTION | eff) --}
{--              (Maybe VHMeas) --}
{-- hSlurpAll [] = do pure Nothing --}
{-- hSlurpAll [p0] = do --}
{--   t0 <- hSlurp p0 --}
{--   let v0 = fst t0 --}
{--   pure v0 --}
{-- hSlurpAll ps = do --}
{--   t0 <- hSlurp $ unsafePartial head ps --}
{--   let v0 = fst t0 --}
{--   vs <- hSlurpAll $ unsafePartial tail ps --}
{--   pure $ v0 -- semigroup not yet implemented... <> vs --}

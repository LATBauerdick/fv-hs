
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Vec where
  import Stuff
  newtype Vec a = Vec { v :: Array Number }

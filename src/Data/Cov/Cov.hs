
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Cov where
  import Stuff
  newtype Cov a = Cov { v :: Array Number }

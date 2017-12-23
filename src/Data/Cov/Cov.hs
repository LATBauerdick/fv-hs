
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Cov where
  import Prelude.Extended
  newtype Cov a = Cov { v :: Array Number }

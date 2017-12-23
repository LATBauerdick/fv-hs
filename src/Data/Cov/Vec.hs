
{-# LANGUAGE DisambiguateRecordFields #-}

module Data.Cov.Vec where
  import Prelude.Extended
  newtype Vec a = Vec { v :: Array Number }

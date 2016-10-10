module Matrix ( inv ) where

import Debug.Trace ( trace )
import Data.Matrix ( inverse, identity, nrows )
import Types ( M )

debug = flip trace

inv :: M -> M
inv m = either invErr id (Data.Matrix.inverse m)
  where invErr s = (identity $ nrows m) `debug` ("🚩" ++ s)



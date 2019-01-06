{-# LANGUAGE RecordWildCards #-}
-- {-# LANGUAGE DeriveFunctor #-}

module Test.Hist ( doHist ) where

import Prelude.Extended

-- import Data.Histogram
-- histo :: (Foldable v, Unbox a, Num a) =>
--          Int
--       -> v Double
--       -> Histogram BinD a
-- histo n v = fillBuilder buildr v
--   where
--     mi = minimum v
--     ma = maximum v
--     bins = binD mi n ma
--     buildr = mkSimple bins

-- doHist = xs where
--   h = histo 4 [1,2,3,5,1,-10,2,3,50,1,6,7,4,6,34,45,20,120,-80]
--   xs = asList


-- import Graphics.Histogram ( plot, histogramNumBins )
-- _ <- plot "cluster-z.png" $ histogramNumBins 90 $ zs vm
-- _ <- plot "cluster-pd.png" $ histogramNumBins 11 $ 1.0 : 0.0 : probs vm

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Backend.Diagrams
import Data.Default.Class
import Data.Colour (opaque)
import Data.Colour.Names (red)
import Control.Lens

chart :: Renderable ()
chart = toRenderable layout
  where
    vals :: [(Double,Double,Double,Double)]
    vals = [ (x,sin (exp x),sin x/2,cos x/10) | x <- [1..20]]
    bars = plot_errbars_values .~ [symErrPoint x y dx dy | (x,y,dx,dy) <- vals]
         $ plot_errbars_title .~"test"
         $ def

    points = plot_points_style .~ filledCircles 2 (opaque red)
           $ plot_points_values .~ [(x,y) |  (x,y,_,_) <- vals]
           $ plot_points_title .~ "test data"
           $ def

    layout = layout_title .~ "Error Bars"
           $ layout_plots .~ [toPlot bars, toPlot points]
           $ def

doHist :: String -> [Number] -> IO ()
doHist s xs = do
  _ <- renderableToFile def{_fo_format=EPS} (s <> ".eps") chart
  pure ()


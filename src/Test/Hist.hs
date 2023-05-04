{-# LANGUAGE RecordWildCards #-}
-- {-# LANGUAGE DeriveFunctor #-}

module Test.Hist ( doHist, doHistXY, doHistVec ) where

import Prelude.Extended hiding ( (.~) )

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
{-
doHistXY :: Text -> [(Number, Number, Number, Number, Number)] -> IO ()
doHistXY s' vals = do
  pure ()
doHist :: Text -> [(Number, Number, Number, Number)] -> IO ()
doHist s' vals = do
  pure ()
doHistVec :: Text -> [(Number, Number, Number, Number, Number)] -> IO ()
doHistVec _ _ = do
  pure ()
-}

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Easy
-- import Graphics.Rendering.Chart.Backend.Cairo
import Graphics.Rendering.Chart.Backend.Diagrams
import Data.Default.Class
import Data.Colour (opaque)
import Data.Colour.Names (red)
import Control.Lens

outProps = fo_size .~ (1024,768)
        $ fo_format .~ SVG
        $ def

pframe = [((x,y),(0.0,0.0)) | x <- range, y <- range] where range = [-0.3, 0.3]
doHistXY :: Text -> [(Number, Number, Number, Number, Number)] -> IO ()
-- plot points with error ellipse
doHistXY s' vals = do
  let s = toString s'
  let chart = toRenderable layout
        where
          bars  = plot_errbars_values .~ [symErrPoint x y dx dy | (x,y,dx,dy,_) <- vals]
                $ plot_errbars_title .~ ( s )
                $ def

          points  = plot_points_style .~ filledCircles 2 (opaque red)
                  $ plot_points_values .~ [(x,y) |  (x,y,_,_,_) <- vals]
                  $ plot_points_title .~ ( s )
                  $ def

          c = opaque blue
          -- vls = vector_line_style . line_color .= c
          -- vhs = vector_head_style . point_color .= c
          frame   = plot_vectors_scale .~ 0.0
                  $ plot_vectors_title .~ ( s )
                  $ plot_vectors_values .~ pframe
                  $ def

          layout = layout_title .~ "FVT"
                 $ layout_x_axis . laxis_generate .~ scaledAxis def (-0.3,0.3)
                 $ layout_y_axis . laxis_generate .~ scaledAxis def (-0.3,0.3)
                 -- $ layout_left_axis_visibility . axis_show_ticks .~ False
                 $ layout_plots .~ [toPlot bars, toPlot points]
                 $ def

  _ <- renderableToFile outProps (s <> ".svg") chart
  pure ()

doHist :: Text -> [(Number, Number, Number, Number)] -> IO ()
-- plot points with error bars
doHist s' vals = do
  let s = toString s'
  let chart = toRenderable layout
        where
          bars  = plot_errbars_values .~ [symErrPoint x y dx dy | (x,y,dx,dy) <- vals]
                $ plot_errbars_title .~ ( s )
                $ def

          points  = plot_points_style .~ filledCircles 2 (opaque red)
                  $ plot_points_values .~ [(x,y) |  (x,y,_,_) <- vals]
                  $ plot_points_title .~ ( s )
                  $ def

          layout = layout_title .~ "FVT"
                 $ layout_plots .~ [toPlot bars, toPlot points]
                 $ def

  _ <- renderableToFile def{_fo_format=SVG} (s <> ".svg") chart
  pure ()

r' x y z = sqrt $ x*x + y*y + z*z
tadd (x1,y1) (x2,y2) = (x1+x2, y1+y2)
efield sign x y = (sign*x/r, sign*y/r) where r = r' x y 10
bfield sign x y = (-sign*y/r/r, sign*x/r/r) where r = r' x y 10
ef (x,y) = efield 1 (x-20) y `tadd` efield (-1) (x+20) y
bf (x,y) = bfield 1 (x-20) y `tadd` bfield (-1) (x+20) y
square a s = [(x,y) | x <- range, y <- range] where range = [-a, -a+s..a] :: [Number]
grid = square 30 3

vvals a s = [((x,y), bf (x,y)) | x <- range, y <- range] where range = [-a, -a+s..a] :: [Number]
vectorField title f grid = fmap plotVectorField $ liftEC $ do
  c <- takeColor
  plot_vectors_mapf .= f
  plot_vectors_grid .= grid
  plot_vectors_style . vector_line_style . line_color .= c
  plot_vectors_style . vector_line_style . line_width .= 0.1
  plot_vectors_style . vector_head_style . point_color .= c
  plot_vectors_title .= title

vectorVal title vals = fmap plotVectorField $ liftEC $ do
  c <- takeColor
  plot_vectors_values .= vals
  plot_vectors_scale .= 0.0
  plot_vectors_style . vector_line_style . line_color .= c
  plot_vectors_style . vector_head_style . point_color .= c
  plot_vectors_title .= title

doHistVec :: Text -> [(Number, Number, Number, Number, Number)] -> IO ()
doHistVec s vals = -- pure ()
  toFile def{_fo_format=SVG} ( (toString s) <> ".svg") $ do
  setColors [ opaque black, opaque blue, opaque red ]

  layout_title .= "Vector Field"
  layout_x_axis . laxis_generate .= scaledAxis def (-0.3,0.3)
  layout_y_axis . laxis_generate .= scaledAxis def (-0.3,0.3)

  -- plot $ vectorPoints "points" vals

  plot $ vectorVal "corr err max" $ fmap ( \(x,y,dx,dy,alf) -> ((x,y),(dx*cos alf,dx*sin alf)) ) vals
  plot $ vectorVal "corr err min" $ fmap ( \(x,y,dx,dy,alf) -> ((x,y),(dy*cos (alf+pi/2.0),dy*sin (alf+pi/2.0))) ) vals
  c <- takeColor
  -- plot_points_style .= filledCircles 2 (opaque red)
  plot $ points "vertices" [(x,y) |  (x,y,_,_,_) <- vals]

  pure ()


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

import Graphics.Rendering.Chart
import Graphics.Rendering.Chart.Easy
import Graphics.Rendering.Chart.Backend.Cairo
import Data.Default.Class
import Data.Colour (opaque)
import Data.Colour.Names (red)
import Control.Lens

outProps = fo_size .~ (1024,768)
        $ fo_format .~ PDF
        $ def

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

          layout = layout_title .~ "FVT"
                 $ layout_plots .~ [toPlot bars, toPlot points]
                 $ def

  _ <- renderableToFile outProps (s <> ".pdf") chart
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

  _ <- renderableToFile def{_fo_format=PDF} (s <> ".pdf") chart
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
  plot_vectors_style . vector_head_style . point_color .= c
  plot_vectors_title .= title

vectorVal title val = fmap plotVectorField $ liftEC $ do
  c <- takeColor
  plot_vectors_values .= val
  plot_vectors_scale .= 0.0
  plot_vectors_style . vector_line_style . line_color .= c
  plot_vectors_style . vector_head_style . point_color .= c
  plot_vectors_title .= title


doHistVec :: Text -> [(Number, Number, Number, Number, Number)] -> IO ()
doHistVec s vals = toFile def{_fo_format=PDF} ( (toString s) <> ".pdf") $ do
  setColors [ opaque blue, opaque black ]

  layout_title .= "Vector Field"
  -- plot $ vectorField "Electric Field" ef grid
  -- plot $ vectorField "Magnetic Field" bf grid
  -- plot $ vectorVal "Mag Field" (vvals 50 5)
  let v = fmap ( \(x,y,dx,dy,alf) -> ((x,y),(dx,dy)) ) vals
  plot $ vectorVal "Mag Field" v

  pure ()


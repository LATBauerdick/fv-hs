{-# LANGUAGE RecordWildCards #-}
-- {-# LANGUAGE DeriveFunctor #-}

module Test.Cluster
  ( doCluster
  , fsmw
  )
where

import           FV.Types                       ( VHMeas(..)
                                                , HMeas(..)
                                                , Prong(..)
                                                , XMeas(..)
                                                , xXMeas
                                                , yXMeas
                                                , XFit(..)
                                                , Chi2(..)
                                                , chi2Vertex
                                                , zVertex
                                                , dzVertex
                                                , rVertex
                                                , drVertex
                                                , z0Helix
                                                )
import           FV.Fit                         ( kAdd
                                                , kAddF
                                                , ksm'
                                                , kChi2
                                                , fit
                                                )
import           Data.Cov.Vec
import           Test.Hist                      ( doHist )
--import Data.Text

import           Prelude.Extended
import           Data.Maybe                     ( mapMaybe
                                                , catMaybes
                                                )
import           Data.List                      (
                                                  sort
                                                -- , foldl'
                                                -- , sortOn
                                                , zip4
                                                , unzip3
                                                )
import qualified Math.Gamma                     ( q )
import           Text.Printf                    ( printf )

doCluster :: VHMeas -> IO Text
doCluster vm' = do
  let vm = cleanup vm'
      v0 = vertex vm' -- beamspot
      s =
        "Test Cluster----------------------------------------------------" <> "\n"
          <> "beam spot -> " <> show v0 <> "\n"
          <> "# tracks  -> " <> show (length <<< helices $ vm') <> "\n"
          <> "# after cleanup-> " <> show (length <<< helices $ vm) <> "\n"
          <> "zs " <> (show . concatMap (unpack <<< to2fix) . zs $ vm) <> "\n"
          <> "rs " <> (show . concatMap (unpack <<< to2fix) . rs $ vm) <> "\n"
      points xs = zip4 (zs xs) (rs xs) (dzs xs) (drs xs)
  _ <- doHist "vertices-z-r" . points $ vm
  return (pack s)

  -- let Node _ ht = vList vm
  -- putStrLn "---------------------------------------------------------"
  -- print $ vList vm
  -- -- print . nProng $ p0
  -- -- print . vertex . measurements $ p0
  -- -- print . fitVertex $ p0
  -- -- print . zip (fitChi2s p0) . map z0Helix . helices . measurements $ p0
  -- case ht of
  --   CEmpty     -> putStrLn "CEmpty"
  --   Node p1 _ -> print $ fitVertex p1

ftvtx :: Prong -> (Number, List Chi2)
ftvtx p = (zVertex . xFit . fitVertex $ p, fitChi2s p)

xFit :: XMeas -> XFit
xFit (XMeas v vv) = XFit v vv (Chi2 1e6)

gamma :: Chi2 -> Number
gamma (Chi2 c2) = Math.Gamma.q 1.0 c2
probs :: VHMeas -> [Number]
probs (VHMeas v hl) =
  filter (> 0.01) $ map (gamma . chi2Vertex . kAddF (xFit v)) hl
zs :: VHMeas -> [Number]
zs (VHMeas v hl) =
  filter (\x -> abs x < 10.0) $ map (zVertex . kAddF (xFit v)) hl
dzs :: VHMeas -> [Number]
dzs (VHMeas v hl) =
  filter (\x -> abs x < 10.0) $ map (dzVertex . kAddF (xFit v)) hl
rs :: VHMeas -> [Number]
rs (VHMeas v hl) =
  filter (\x -> abs x < 10.0) $ map (rVertex . kAddF (xFit v)) hl
drs :: VHMeas -> [Number]
drs (VHMeas v hl) =
  filter (\x -> abs x < 10.0) $ map (drVertex . kAddF (xFit v)) hl

cleanup :: VHMeas -> VHMeas
-- remove vm helices that are incompatible with vm vertex
cleanup (VHMeas v hl) = VHMeas v hl'
  where hl' = sortOn z0Helix . mapMaybe (filtProb 0.01 v) $ hl

filtProb :: Number -> XMeas -> HMeas -> Maybe HMeas
filtProb cut v h = mh
 where
-- check chi2 of this helix w/r to vertex position v
  vf         = kAddF (xFit v) h
  zvf        = zVertex vf
  Chi2 chi2f = chi2Vertex vf
  prob       = gamma $ Chi2 (chi2f / 2.0) -- chi2 distribution with NDOF=2
  good       = (prob > cut) && abs zvf < 10.0
  mh         = if good then Just h else Nothing

-- Now the "v" parameter can be mapped over without any care for list invariants
data HList v = CEmpty | Node v (HList v) deriving (Show, Functor)

vList :: VHMeas -> HList Prong
vList vm = Node p vRight
 where
  (p, vmr) = cluster vm
  vRight   = case vmr of
    Nothing  -> CEmpty
    Just vm' -> vList vm'

wght :: Number -> Chi2 -> Double -- weight function with Temperature t
wght t (Chi2 chi2) = w
 where
  chi2cut = 3.0
  w       = 1.0 / (1.0 + exp ((chi2 - chi2cut) / 2.0 / t))

-- cluster'' :: VHMeas -> (Prong, Maybe VHMeas)
-- cluster'' vm | trace ("--> cluster called with " <> (show . length . helices $ vm) <> " helices, initial vertex at " <> (show . vertex $ vm) ) False = undefined
-- cluster'' (VHMeas v hl) = trace (
--         "--> cluster debug:"
--         <> "\n--> cluster fit zs " <> show zs
--         <> "\n--> cluster fit vertex " <> show (fitVertex p)
--         <> "\n--> cluster # remaining helices " <> show ( length hlr) <> (show . fmap z0Helix . take 3 $ hlr)
--         ) ( p, r ) where
--   -- split off 10 h with lowest z and determine initial fit position with FSMW
--   (hl', hlr) = splitAt 10 . sortOn z0Helix $ hl

--   -- constract vertex v0 at z=z0 as starting point, using cov matrix from v
--   zs = sort <<< filter (\x -> abs x < 10.0) <<< map (zVertex <<< kAddF (xFit v)) $ hl'
--   z0 = fsmw (length zs) zs
--   XMeas _ cx0  = v
--   v0 = XMeas Vec {v= fromList [xXMeas v, yXMeas v, z0]} cx0

--   -- do the fit with v0
--   p             = fit (VHMeas v0 hl')

--   r = case length hlr of
--         0 -> Nothing
--         _ -> Just (VHMeas v hlr)

cluster :: VHMeas -> (Prong, Maybe VHMeas)
cluster vm
  | trace
    (  "--> cluster called with "
    <> (tshow . length . helices $ vm)
    <> " helices, initial vertex at "
    <> (tshow . vertex $ vm)
    )
    False
  = undefined
cluster (VHMeas v hllll) = trace (
  pack (  "--> cluster debug:"
  <> "\n--> cluster zs=" <> take 160 (show zs)
  <> "\n--> cluster z0=" <> unpack (to3fix (z0 :: Number))
  <> "\n--> cluster v0=" <> show v0
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 10 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 20 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 30 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 40 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 50 hl)
  <> "\n--> cluster fit v0 hl0 " <> (show . ftvtx . fit $ VHMeas v0 $ take 10 $ drop 60 hl)
  <> "\n--> cluster v1: #hl0=" <> (show . length . catMaybes $ hl0) <> " v1=" <> show v1
  <> "\n--> cluster chi2s1" <> (take 160 . foldl (\s c -> s <> bDouble c) "") c21s
  <> "\n--> cluster weights1" <> (take 160 . foldl (\s w -> s <> sDouble w) "") ws1
  <> "\n--> cluster v2=" <> show v2
  <> "\n--> cluster chi2s2" <> (take 160 . foldl (\s c -> s <> bDouble c) "") c22s
  <> "\n--> cluster weights2" <> (take 160 . foldl (\s w -> s <> sDouble w) "") ws2
  <> "\n--> cluster v3=" <> show v3
  <> "\n--> cluster chi2s3" <> (take 160 . foldl (\s c -> s <> bDouble c) "") c23s
  <> "\n--> cluster weights3" <> (take 160 . foldl (\s w -> s <> sDouble w) "") ws3
  <> "\n--> cluster v4=" <> show v4
  <> "\n--> cluster chi2s4" <> (take 160 . foldl (\s c -> s <> bDouble c) "") c24s
  <> "\n--> cluster weights4" <> (take 160 . foldl (\s w -> s <> sDouble w) "") ws4
  <> "\n--> cluster v5=" <> show v5
  <> "\n-----------------------------------------------------------------"
  <> "\n--> cluster #hl1=" <> (show . length . catMaybes $ hl1) <> " vv=" <> show vv
  <> printf "\n--> cluster nothing=%5d just=%5d" (length hlnothing) (length hljust)
  <> "\n--> cluster nProng=" <> show nProng
  ))
  (p, r)
 where

  sDouble :: Double -> [Char]
  sDouble c = if c > 0.001 then unpack (to3fix c) else " eps"

  bDouble :: Chi2 -> [Char]
  bDouble (Chi2 c) = if c < 999.9 then unpack (to1fix c) else " big"

  hl          = hllll

  -- FSMW method to find initial z for primary vertex
  zs = sort . filter (\x -> abs x < 10.0) . map (zVertex . kAddF (xFit v)) $ hl
  z0          = fsmw (length zs) zs
  -- constract vertex v0 at z=z0 as starting point, using cov matrix from v
  XMeas _ cx0 = v
  v0          = XMeas Vec {v = fromList [xXMeas v, yXMeas v, z0]} cx0
  -- filter hl for v1, with starting point v0
  hl0 :: [Maybe HMeas]
  hl0 = map (filtProb 0.00001 v0) hl
  v1  = foldl
    (\v_ mh -> case mh of
      Just h  -> kAdd v_ h
      Nothing -> v_
    )
    v0
    hl0

  kAddW' :: XMeas -> (Maybe HMeas, Number) -> XMeas
  kAddW' v_ (mh, _) = case mh of
    Just h  -> kAdd v_ h
    Nothing -> v_
  annealingSchedule :: NonEmpty Number;
  annealingSchedule = 256.0 :| [64.0, 16.0, 4.0, 1.0]
  t0                = head annealingSchedule
  -- t1 = annealingSchedule !! 3

  c21s              = map
    (\mh -> case mh of
      Just h  -> kChi2 v1 h
      Nothing -> Chi2 0.0
    )
    hl0
  ws1  = fmap (max 0.001 . wght t0) c21s
  v2   = foldl kAddW' v1 $ zip hl0 ws1
  c22s = map
    (\mh -> case mh of
      Just h  -> kChi2 v2 h
      Nothing -> Chi2 0.0
    )
    hl0
  ws2  = fmap (max 0.001 . wght t0) c22s
  v3   = foldl kAddW' v2 $ zip hl0 ws2
  c23s = map
    (\mh -> case mh of
      Just h  -> kChi2 v3 h
      Nothing -> Chi2 0.0
    )
    hl0
  ws3  = fmap (max 0.001 . wght 4.0) c23s
  v4   = foldl kAddW' v3 $ zip hl0 ws3
  c24s = map
    (\mh -> case mh of
      Just h  -> kChi2 v4 h
      Nothing -> Chi2 0.0
    )
    hl0
  ws4 = fmap (max 0.001 . wght 1.0) c24s
  v5  = foldl kAddW' v4 $ zip hl0 ws4

  vv0 = v5

  -- cut any track with prob < 0.001 w/r to v1
  hl1 = map (filtProb 0.00001 vv0) hl
  -- re-do filter with initial position and filtered helices
  vv  = foldl
    (\v_ h -> case h of
      Just h' -> kAdd v_ h'
      Nothing -> v_
    )
    vv0
    hl1
  -- smooth with vv
  ll = zip (map (ksm' vv) hl1) hl
  hlnothing :: [HMeas]
  hlnothing = [ h | (Nothing, h) <- ll ]
  hljust :: [HMeas]
  (qljust, c2just, hljust) = unzip3 [ (q, c, h) | (Just (q, c), h) <- ll ]

  fitVertex                = vv
  fitMomenta               = qljust
  fitChi2s                 = c2just
  nProng                   = length qljust
  measurements             = VHMeas v hljust
  p                        = Prong {..}

  nnn                      = fromIntegral (length hlnothing) `div` 3
  r = if null hlnothing then Nothing else Just (VHMeas v (drop nnn hlnothing))

-- -------------------------------------------------------
-- use robust Mode finding to get initial vertex position
-- -------------------------------------------------------
-- Fraction-of Sample Mode with Weight method,
-- see Frühwirth & Waltenberger CMS Note 2007/008
--newtype WeightedPoint = WeightedPoint Number
type WeightedPoint = Number
fsmw :: Int -> List WeightedPoint -> WeightedPoint
fsmw 0 []            = error $ "Test.Cluster.fsmw got empty list"
fsmw 1 [a]           = a -- head xs
fsmw 2 [x0, x1]     = 0.5 * (x1 + x0)
fsmw 3 [x0, x1, x2] = case 2 * x1 - x0 - x2 of
  xx | xx < 0 -> (x0 + x1) / 2.0
  xx | xx > 0 -> (x1 + x2) / 2.0
  _           -> x1
-- fsmw 4 [x0, x1, x2, x3] = x where
--   w  = weightOfInterval [x0, x1, x2, x3]
--   w0 = weightOfInterval [x0, x1, x2]
--   w1 = weightOfInterval [x1, x2, x3]
--   x  = if w0/w < w1/w then fsmw 3 [x0, x1, x2] else fsmw 3 [x1, x2, x3] `debug` (show w <> ", " <> show w0 <> ", " <> show w1)
fsmw n xs = h
 where
  h = if n > 6
    then hsm n xs -- for large n, fall back to no-weight hsm formula
    else fsmw n' xs'
   where
    alpha = 0.5 :: Number
    n'    = ceiling (fromIntegral n * alpha)
    findMin
      :: (WeightedPoint, Int) -> (WeightedPoint, Int) -> (WeightedPoint, Int)
    findMin (w0, j0) (w, j) = if w < w0 || w0 < 0 then (w, j) else (w0, j0)
    (_, j') =
      foldl' findMin (-1.0, 0)
        . zip
            (map (\j -> weightOfInterval . take n' . drop j $ xs) [0 .. n - n'])
        $ [0 ..]
    xs' = take n' . drop j' $ xs -- `debug` ("fsmw--> " <> show n <> ", " <> show n' <> ", " <> show j' <> ", " <> show xs)

-- HalfSampleMode
-- see Bickel & Frühwirth, On a Fast, Robust Estimator of the Mode, 2006
-- http://ideas.repec.org/a/eee/csdana/v50y2006i12p3500-3530.html
hsm :: Int -> List WeightedPoint -> WeightedPoint
hsm n xs = fsmw n' xs'
 where
  alpha = 0.5 :: Number
  n'    = ceiling (fromIntegral n * alpha)
  xns   = drop (n' - 1) xs
  wmin  = unsafeLast xs - unsafeHead xs
  findMin
    :: (WeightedPoint, Int) -> (WeightedPoint, Int) -> (WeightedPoint, Int)
  findMin (w0, j0) (w, j) = if w < w0 then (w, j) else (w0, j0)
  calcWeight :: WeightedPoint -> WeightedPoint -> WeightedPoint
  calcWeight = (-)
  (_, j') = foldl' findMin (wmin, 0) . zip (zipWith calcWeight xns xs) $ [0 ..]
  xs' = take n' . drop j' $ xs --`debug` ("hsm---> " <> show n <> ", " <> show n' <> ", " <> show j' <> ", " <> show xs <> ", " <> show xns)

wDist :: WeightedPoint -> WeightedPoint -> WeightedPoint
wDist x0 x1 = 1.0 / sqrt (x1 - x0 + dmin) where dmin = 0.001 -- 10 µm
weightOfInterval :: List WeightedPoint -> WeightedPoint
weightOfInterval xs = w
 where --`debug` ("weightOfInterval--> " <> show xs <> ", " <> show ws <> ", " <> show w) where
  ws = [ wDist x0 x1 | x0 <- xs, x1 <- xs, x1 > x0 ]
  w  = (unsafeLast xs - unsafeHead xs) / sum ws




{-
-- gamma function P(a, x) from https://wiki.haskell.org/Gamma_and_Beta_function
-- approximation is taken from [Lanczos, C. 1964 SIAM Journal on Numerical Analysis,
-- ser. B, vol. 1, pp. 86-96]

-- standard normal CDF https://www.johndcook.com/blog/haskell-erf/
erf :: Number -> Number
erf x = sign*y
    where
        a1 =  0.254829592
        a2 = -0.284496736
        a3 =  1.421413741
        a4 = -1.453152027
        a5 =  1.061405429
        p  =  0.3275911

        -- Abramowitz and Stegun formula 7.1.26
        sign = if x > 0
                   then  1
                   else -1
        t  =  1.0/(1.0 + p* abs x)
        y  =  1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)

test_erf :: Bool
test_erf = maximum [ abs(erf x - y)  | (x, y) <- zip xs ys ] < epsilon
    where
        epsilon = 1.5e-7 -- accuracy promised by A&S
        xs = [-3, -1, 0.0, 0.5, 2.1 ]
        ys = [-0.999977909503,
              -0.842700792950,
               0.0,
               0.520499877813,
               0.997020533344]

-- standard normal CDF https://www.johndcook.com/blog/haskell-phi/
phi :: Number -> Number
phi x = y
    where
        a1 =  0.254829592
        a2 = -0.284496736
        a3 =  1.421413741
        a4 = -1.453152027
        a5 =  1.061405429
        p  =  0.3275911

        -- Abramowitz and Stegun formula 7.1.26
        sign = if x > 0
                   then  1
                   else -1
        t = 1.0/(1.0 + p * abs x / sqrt 2.0)
        e = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x/2.0)
        y = 0.5*(sign*e + 1.0)

test_phi :: Bool
test_phi = maximum [ abs(phi x - y) | (x, y) <- zip xs ys ] < epsilon
    where
        epsilon = 1.5e-7 -- accuracy promised by A&S
        xs = [-3, -1, 0.0, 0.5, 2.1 ]
        ys = [0.00134989803163,
              0.158655253931,
              0.5,
              0.691462461274,
              0.982135579437]


-- returns the number of samples and the mean,
-- but does not directly return variance, skewness and kurtosis.
-- Instead it returns moments from which these statistics can easily be calculated
-- using the mvks function
moments (n, m1, m2, m3, m4) x = (n', m1', m2', m3', m4')
        where
            n' = n + 1
            delta = x - m1
            delta_n = delta / n'
            delta_n2 = delta_n**2
            t = delta*delta_n*n
            m1' = m1 + delta_n
            m4' = m4 + t*delta_n2*(n'*n' - 3*n' + 3) + 6*delta_n2*m2 - 4*delta_n*m3
            m3' = m3 + t*delta_n*(n' - 2) - 3*delta_n*m2
            m2' = m2 + t

mvsk (n, m1, m2, m3, m4) = (m1, m2/(n-1.0), (sqrt n)*m3/m2**1.5, n*m4/m2**2 - 3.0)
-}

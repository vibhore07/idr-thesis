/**** 
NEIGHBOURHOOD PLANNING ONLY (Indore Metro station areas)
Creates thesis-ready neighbourhood planning analysis/proposal maps around given stations.

What this script outputs (Drive exports):
NBH_1_WalkCatchments: 5/10/15-minute walk catchments (cost-distance proxy)
NBH_2_GreenBlueNetwork: green + water + suggested green connectors within walk10
NBH_3_BuiltIntensity: built-up density heatmap within walk10 (proxy for intensity/pressure)
NBH_4_HeatComfort: LST heat proxy + green/water + walk10 overlay (thermal comfort planning)

Compatible with your runtime:
- Uses instance method: costImage.cumulativeCost(source, maxDistance, geodeticDistance)
- No named args.

IMPORTANT LIMITATION:
- True network-based isochrones require OSM road network analysis outside GEE.
  This is a robust GEE-only approximation suitable for thesis mapping + proposals.

Date: 2026-05-02
****/

// -------------------------
// 0) Metro stations (given)
// -------------------------
var stations = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([75.797608, 22.739112]), {name: 'Devi Ahilya Bai Holkar Terminal'}),
  ee.Feature(ee.Geometry.Point([75.801004, 22.745436]), {name: 'Maharani Lakshmi Bai'}),
  ee.Feature(ee.Geometry.Point([75.804299, 22.754001]), {name: 'Rani Avanti Bai Lodhi'}),
  ee.Feature(ee.Geometry.Point([75.809859, 22.761346]), {name: 'Super Corridor (Barawagrada)'}),
  ee.Feature(ee.Geometry.Point([75.819195, 22.771952]), {name: 'Veerangana Jhalkari Bai (SC-03)'}),
  ee.Feature(ee.Geometry.Point([75.826110, 22.778380]), {name: 'Super Corridor 2'}),
  ee.Feature(ee.Geometry.Point([75.837161, 22.786602]), {name: 'Super Corridor 1'}),
  ee.Feature(ee.Geometry.Point([75.846651, 22.789829]), {name: 'Bhawarsala Square'})
]);

// Study area (neighbourhood planning context)
var padMeters = 6000;
var studyArea = stations.geometry().buffer(padMeters).bounds();
Map.centerObject(studyArea, 13);

// Neighbourhood focus radius for proposals (10-min walk)
var thr5  = 450;   // ~5 min effective meters
var thr10 = 900;   // ~10 min effective meters
var thr15 = 1350;  // ~15 min effective meters

// cumulativeCost config
var maxDistance = 20000;
var geodeticDistance = true;

// -------------------------
// 1) Base imagery (Sentinel-2)
// -------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2025-10-01', '2026-03-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .median()
  .clip(studyArea);

var s2Vis = {bands: ['B4','B3','B2'], min: 0, max: 3000, gamma: 1.15};

// -------------------------
// 2) Land cover (WorldCover v200)
// -------------------------
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);

// Masks
var built = wc.eq(50).selfMask();
var green = wc.eq(10).or(wc.eq(20)).or(wc.eq(30)).or(wc.eq(40)).selfMask();
var water = wc.eq(80).selfMask();

// Built-up density (proxy for neighbourhood intensity)
var builtDensity = built.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.circle({radius: 500, units: 'meters', normalize: true})
}).rename('built_density');

// -------------------------
// 3) Walking friction surface (cost image)
// -------------------------
var dem = ee.Image('USGS/SRTMGL1_003').clip(studyArea);
var slope = ee.Terrain.slope(dem);

// Landcover friction (tunable)
var frictionLC = ee.Image(1.0)
  .where(wc.eq(50), 1.2)     // built-up
  .where(wc.eq(60), 1.1)     // barren
  .where(wc.eq(40), 1.0)     // cropland
  .where(wc.eq(10), 1.0)     // trees
  .where(wc.eq(30), 1.0)     // grass
  .where(wc.eq(80), 100.0);  // water ~ barrier

// Slope penalty
var frictionSlope = slope.multiply(0.03).add(1.0);

// Final cost surface
var friction = frictionLC.multiply(frictionSlope).rename('friction');

// -------------------------
// 4) Walk catchments (Neighbourhood structure)
// -------------------------
var sourceImg = ee.Image().byte().paint(stations, 1); // sources at stations

// Correct for your runtime: costImage.cumulativeCost(source, maxDistance, geodeticDistance)
var costDist = friction.cumulativeCost(sourceImg, maxDistance, geodeticDistance).rename('costDist');

var walk5  = costDist.lte(thr5).selfMask();
var walk10 = costDist.lte(thr10).selfMask();
var walk15 = costDist.lte(thr15).selfMask();

// -------------------------
// 5) Heat proxy (Landsat 9 LST) for comfort planning
// -------------------------
function maskL9(img) {
  var qa = img.select('QA_PIXEL');
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  return img.updateMask(cloud.not()).updateMask(shadow.not());
}

var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2025-03-01', '2025-06-30') // warm season window
  .map(maskL9)
  .median()
  .clip(studyArea);

// LST Kelvin -> Celsius
var lstC = l9.select('ST_B10')
  .multiply(0.00341802).add(149.0)
  .subtract(273.15)
  .rename('LST_C');

var lstSmooth = lstC.focal_mean({radius: 300, units: 'meters'}).rename('LST_C_smooth');

// -------------------------
// 6) Proposal logic layers (neighbourhood planning)
// -------------------------

// A) Green-blue access: green+water within 10-min
var greenBlue = green.unmask(0).add(water.unmask(0)).gt(0).selfMask();

// B) Suggested green connectors (simple "where to connect" proxy):
// We compute distance-to-greenBlue; areas far from green are candidates for new street trees/pocket parks.
// (This is a planning heuristic map.)
var distToGreen = greenBlue
  .fastDistanceTransform(30).sqrt()  // in pixels
  .multiply(10)                      // approx meters (WorldCover ~10m)
  .rename('dist_to_green_m');

// High-need zones: > 300m from green/blue within walk10
var greenNeed = distToGreen.gt(300).selfMask().updateMask(walk10);

// -------------------------
// 7) Styling
// -------------------------
var stationsVis = stations.style({color: 'ffffff', pointSize: 6, width: 2, fillColor: 'd50000'});
var boundary = ee.Image().paint(studyArea, 1, 2).visualize({palette: ['000000'], opacity: 1});

var walk5Vis  = walk5.visualize({palette: ['00c853'], opacity: 0.18});
var walk10Vis = walk10.visualize({palette: ['00c853'], opacity: 0.35});
var walk15Vis = walk15.visualize({palette: ['00c853'], opacity: 0.12});

var greenVis = green.visualize({palette: ['1b5e20'], opacity: 0.65});
var waterVis = water.visualize({palette: ['1565c0'], opacity: 0.80});

var builtDenVis = builtDensity.visualize({
  min: 0, max: 0.6,
  palette: ['ffffff','ffe082','ffb300','fb8c00','e65100','3e2723']
});

var heatVis = lstSmooth.visualize({
  min: 25, max: 45,
  palette: ['2c7bb6','abd9e9','ffffbf','fdae61','d7191c']
});

var greenNeedVis = greenNeed.visualize({palette: ['ff1744'], opacity: 0.40}); // red = green deficit

// -------------------------
// 8) Final neighbourhood planning maps (export images)
// -------------------------

// NBH 1: Walk catchments (neighbourhood structure)
var NBH_1_WalkCatchments = s2.visualize(s2Vis)
  .blend(walk15Vis).blend(walk10Vis).blend(walk5Vis)
  .blend(stationsVis)
  .blend(boundary)
  .clip(studyArea);

// NBH 2: Green-Blue network + deficit zones within walk10 (proposal guidance)
var NBH_2_GreenBlueNetwork = s2.visualize(s2Vis)
  .blend(walk10Vis)
  .blend(greenVis).blend(waterVis)
  .blend(greenNeedVis)       // where to propose pocket parks / shaded streets
  .blend(stationsVis)
  .blend(boundary)
  .clip(studyArea);

// NBH 3: Built intensity within walk10 (pressure / redevelopment focus)
var walk10Mask = ee.Image().byte().paint(walk10.geometry(), 1).selfMask();
var NBH_3_BuiltIntensity = builtDenVis
  .updateMask(walk10Mask)
  .unmask(s2.visualize(s2Vis))
  .blend(walk10Vis)
  .blend(stationsVis)
  .blend(boundary)
  .clip(studyArea);

// NBH 4: Heat comfort + green/water + walk10
var NBH_4_HeatComfort = heatVis
  .blend(greenVis).blend(waterVis)
  .blend(walk10Vis)
  .blend(stationsVis)
  .blend(boundary)
  .clip(studyArea);

// -------------------------
// 9) Preview in GEE
// -------------------------
Map.addLayer(NBH_1_WalkCatchments, {}, 'EXPORT: NBH_1 Walk Catchments', true);
Map.addLayer(NBH_2_GreenBlueNetwork, {}, 'EXPORT: NBH_2 Green-Blue + Deficit', false);
Map.addLayer(NBH_3_BuiltIntensity, {}, 'EXPORT: NBH_3 Built Intensity (Walk10)', false);
Map.addLayer(NBH_4_HeatComfort, {}, 'EXPORT: NBH_4 Heat Comfort', false);

// Debug (optional)
Map.addLayer(costDist, {min: 0, max: 3000, palette: ['ffffff','ccebc5','7bccc4','2b8cbe','08589e']}, 'Debug: costDist', false);
Map.addLayer(distToGreen, {min: 0, max: 1000, palette: ['1b5e20','ffffbf','ff1744']}, 'Debug: distToGreen(m)', false);

// -------------------------
// 10) Exports to Google Drive
// -------------------------
var folder = 'GEE_Exports';
var scaleS2 = 10;
var scaleL9 = 30;

Export.image.toDrive({
  image: NBH_1_WalkCatchments,
  description: 'Indore_NBH_1_WalkCatchments',
  folder: folder,
  fileNamePrefix: 'indore_nbh_1_walk_catchments',
  region: studyArea,
  scale: scaleS2,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: NBH_2_GreenBlueNetwork,
  description: 'Indore_NBH_2_GreenBlueNetwork',
  folder: folder,
  fileNamePrefix: 'indore_nbh_2_green_blue_network',
  region: studyArea,
  scale: scaleS2,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: NBH_3_BuiltIntensity,
  description: 'Indore_NBH_3_BuiltIntensity',
  folder: folder,
  fileNamePrefix: 'indore_nbh_3_built_intensity_walk10',
  region: studyArea,
  scale: scaleS2,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: NBH_4_HeatComfort,
  description: 'Indore_NBH_4_HeatComfort',
  folder: folder,
  fileNamePrefix: 'indore_nbh_4_heat_comfort',
  region: studyArea,
  scale: scaleL9,
  maxPixels: 1e13
});
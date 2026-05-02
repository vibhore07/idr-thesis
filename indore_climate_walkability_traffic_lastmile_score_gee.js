/****
Indore Metro (given stations): Climate + Walkability + Traffic-risk + Last-mile connectivity score
Site: same station set (Super Corridor corridor zone)

What this script produces (export GeoTIFFs):
C1_UHI_HeatRisk: Land Surface Temperature (Landsat 9) heat risk + station zones
C2_GreenShadeIndex: NDVI (Sentinel-2) as shade/green proxy + deficit zones
C3_TrafficRiskProxy: Nighttime lights (VIIRS) + built density as a traffic/activity proxy
C4_WalkabilityScore: Composite walkability suitability (green, heat, slope, built pressure)
C5_LastMileConnectivityScore: Composite last-mile score per pixel (proximity to stations + walkability + activity)
C6_PriorityInterventions: classes for where to prioritize (shade, crossings, traffic calming, feeders)

Important (honest limitation):
- GEE cannot compute true traffic volumes/speeds without local traffic count data.
- We use robust proxies: VIIRS night lights + built-up density + proximity to stations + slope + heat + greenery.
- If you later provide actual traffic count points/road network, we can replace the proxy layer.

NO cumulativeCost used (stable).
Date: 2026-05-02
****/

// -------------------------
// 0) Stations (given)
// -------------------------
var stations = ee.FeatureCollection([
  ee.Feature(ee.Geometry.Point([75.797608, 22.739112]), {name: 'Devi Ahilya Bai Holkar Terminal', tier: 'primary'}),
  ee.Feature(ee.Geometry.Point([75.801004, 22.745436]), {name: 'Maharani Lakshmi Bai', tier: 'secondary'}),
  ee.Feature(ee.Geometry.Point([75.804299, 22.754001]), {name: 'Rani Avanti Bai Lodhi', tier: 'secondary'}),
  ee.Feature(ee.Geometry.Point([75.809859, 22.761346]), {name: 'Super Corridor (Barawagrada)', tier: 'primary'}),
  ee.Feature(ee.Geometry.Point([75.819195, 22.771952]), {name: 'Veerangana Jhalkari Bai (SC-03)', tier: 'secondary'}),
  ee.Feature(ee.Geometry.Point([75.826110, 22.778380]), {name: 'Super Corridor 2', tier: 'secondary'}),
  ee.Feature(ee.Geometry.Point([75.837161, 22.786602]), {name: 'Super Corridor 1', tier: 'secondary'}),
  ee.Feature(ee.Geometry.Point([75.846651, 22.789829]), {name: 'Bhawarsala Square', tier: 'primary'})
]);

// -------------------------
// 1) Study area + station influence rings
// -------------------------
var r_core = 400, r_trans = 800, r_infl = 1500;
var pad = 9000;
var studyArea = stations.geometry().buffer(pad).bounds();
Map.centerObject(studyArea, 13);

// Union influence for masking outputs
var inflUnion = stations.map(function(f){ return f.buffer(r_infl); }).geometry().dissolve();
var mask_infl = ee.Image().byte().paint(inflUnion, 1).selfMask();

// Station proximity (distance raster)
var stationSeeds = ee.Image().byte().paint(stations, 1).selfMask();
var distToStation = stationSeeds.fastDistanceTransform(30).sqrt().multiply(10).rename('dist_to_station_m');

// Normalize proximity score (1 near, 0 far at 1500m)
var prox = ee.Image(1).subtract(distToStation.unitScale(0, r_infl).clamp(0, 1)).rename('prox_score').updateMask(mask_infl);

// -------------------------
// 2) Climate layers: Heat (LST) + NDVI (green/shade proxy)
// -------------------------

// (A) Landsat 9 LST (heat risk proxy)
function maskL9(img) {
  var qa = img.select('QA_PIXEL');
  var cloud = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  return img.updateMask(cloud.not()).updateMask(shadow.not());
}

var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(studyArea)
  .filterDate('2025-03-01', '2025-06-30')
  .map(maskL9)
  .median()
  .clip(studyArea);

var lstC = l9.select('ST_B10')
  .multiply(0.00341802).add(149.0)
  .subtract(273.15)
  .rename('LST_C');

var lstSmooth = lstC.focal_mean({radius: 300, units: 'meters'}).rename('LST_C_smooth');

// Heat risk score (0 cool .. 1 hot)
var heatScore = lstSmooth.unitScale(25, 45).clamp(0, 1).rename('heat_score').updateMask(mask_infl);

// (B) Sentinel-2 NDVI (shade/green proxy)
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2025-10-01', '2026-03-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .median()
  .clip(studyArea);

var s2Vis = {bands: ['B4','B3','B2'], min: 0, max: 3000, gamma: 1.15};

var ndvi = s2.normalizedDifference(['B8','B4']).rename('NDVI');
var greenScore = ndvi.unitScale(0.15, 0.65).clamp(0, 1).rename('green_score').updateMask(mask_infl); // higher = greener

// Green deficit (for proposals)
var greenDeficit = greenScore.lte(0.35).selfMask().updateMask(mask_infl);

// -------------------------
// 3) Terrain barrier: slope penalty (walkability friction)
// -------------------------
var dem = ee.Image('USGS/SRTMGL1_003').clip(studyArea);
var slope = ee.Terrain.slope(dem).rename('slope_deg');

// Slope score: flatter = better (0 steep .. 1 flat)
var slopeScore = ee.Image(1).subtract(slope.unitScale(0, 10).clamp(0, 1)).rename('slope_score').updateMask(mask_infl);

// -------------------------
// 4) Built / activity proxies (walk demand + traffic risk proxy)
// -------------------------

// (A) Built-up density from WorldCover
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);
var built = wc.eq(50).selfMask();
var builtDensity = built.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.circle({radius: 500, units: 'meters', normalize: true})
}).rename('built_density').updateMask(mask_infl);

// Built pressure score (0 low .. 1 high)
var builtScore = builtDensity.unitScale(0, 0.6).clamp(0, 1).rename('built_score').updateMask(mask_infl);

// (B) Nighttime lights (VIIRS) as traffic/activity intensity proxy
// Use monthly composites to reduce noise
var viirs = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
  .filterBounds(studyArea)
  .filterDate('2025-01-01', '2025-12-31')
  .median()
  .select('avg_rad')
  .clip(studyArea);

var lightsScore = viirs.unitScale(0, 30).clamp(0, 1).rename('lights_score').updateMask(mask_infl);

// Traffic risk proxy: combine lights + built pressure
var trafficRisk = lightsScore.multiply(0.6).add(builtScore.multiply(0.4))
  .rename('traffic_risk_proxy')
  .updateMask(mask_infl);

// -------------------------
// 5) WALKABILITY SCORE (0..100)
// -------------------------
// Logic: high green + flat + near station + low heat + manageable traffic risk
// (Traffic risk is treated as negative for pedestrian comfort/safety)

var walkability01 = greenScore.multiply(0.30)
  .add(slopeScore.multiply(0.15))
  .add(prox.multiply(0.20))
  .add(ee.Image(1).subtract(heatScore).multiply(0.20))
  .add(ee.Image(1).subtract(trafficRisk).multiply(0.15))
  .rename('walkability_01')
  .updateMask(mask_infl);

var walkability100 = walkability01.multiply(100).rename('walkability_score_0_100');

// -------------------------
// 6) LAST-MILE CONNECTIVITY SCORE (0..100)
// -------------------------
// Definition (thesis-ready):
// Last-mile score = proximity to station + walkability + activity (built) - heat penalty - traffic risk penalty

var lastMile01 = prox.multiply(0.35)
  .add(walkability01.multiply(0.35))
  .add(builtScore.multiply(0.15))       // demand + viability for feeders/NMT
  .add(ee.Image(1).subtract(heatScore).multiply(0.10))
  .add(ee.Image(1).subtract(trafficRisk).multiply(0.05))
  .rename('lastmile_01')
  .updateMask(mask_infl);

var lastMile100 = lastMile01.multiply(100).rename('lastmile_score_0_100');

// -------------------------
// 7) PRIORITY INTERVENTION CLASSES (proposal map)
// -------------------------
// Class meanings:
// 1 = Shade + trees + cooling priority (hot + green deficit)
// 2 = Traffic calming + safe crossings priority (high traffic risk + low walkability)
// 3 = Last-mile feeder/NMT network priority (high demand + low prox / low last-mile)
// 4 = Maintain/upgrade (already good score)

var hotGreenDef = heatScore.gte(0.7).and(greenScore.lte(0.35));
var highTrafficLowWalk = trafficRisk.gte(0.6).and(walkability01.lte(0.45));
var demandLowAccess = builtScore.gte(0.35).and(lastMile01.lte(0.45));

var interventionClass = ee.Image(4)
  .where(hotGreenDef, 1)
  .where(highTrafficLowWalk, 2)
  .where(demandLowAccess, 3)
  .updateMask(mask_infl)
  .rename('intervention_class');

// -------------------------
// 8) Map compositions (export-ready visuals)
// -------------------------
var boundaryImg = ee.Image().paint(studyArea, 1, 2).visualize({palette: ['000000'], opacity: 1});
var stationsImg = stations.style({color: 'ffffff', pointSize: 6, width: 2, fillColor: 'd50000'});

// C1 Heat risk
var C1_UHI_HeatRisk = s2.visualize(s2Vis)
  .blend(lstSmooth.visualize({min: 25, max: 45, palette: ['2c7bb6','abd9e9','ffffbf','fdae61','d7191c'], opacity: 0.75}))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// C2 Green/shade index
var C2_GreenShadeIndex = s2.visualize(s2Vis)
  .blend(ndvi.visualize({min: 0.0, max: 0.7, palette: ['8d6e63','fff59d','66bb6a','1b5e20'], opacity: 0.75}))
  .blend(greenDeficit.visualize({palette: ['ff1744'], opacity: 0.25}))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// C3 Traffic risk proxy
var C3_TrafficRiskProxy = s2.visualize(s2Vis)
  .blend(trafficRisk.visualize({min: 0, max: 1, palette: ['1b5e20','ffffbf','ff6f00','b71c1c'], opacity: 0.75}))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// C4 Walkability score
var C4_WalkabilityScore = s2.visualize(s2Vis)
  .blend(walkability100.visualize({min: 0, max: 100, palette: ['b71c1c','ffb300','fff59d','66bb6a','1b5e20'], opacity: 0.80}))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// C5 Last-mile connectivity score
var C5_LastMileConnectivityScore = s2.visualize(s2Vis)
  .blend(lastMile100.visualize({min: 0, max: 100, palette: ['b71c1c','ffb300','fff59d','66bb6a','1b5e20'], opacity: 0.80}))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// C6 Priority interventions
var C6_PriorityInterventions = s2.visualize(s2Vis)
  .blend(interventionClass.visualize({
    min: 1, max: 4,
    palette: [
      '1e88e5', // 1 shade/cooling (blue)
      'e53935', // 2 traffic calming (red)
      'fb8c00', // 3 feeder/NMT (orange)
      '43a047'  // 4 maintain (green)
    ],
    opacity: 0.75
  }))
  .blend(stationsImg)
  .blend(boundaryImg)
  .clip(studyArea);

// -------------------------
// 9) Preview layers
// -------------------------
Map.addLayer(C4_WalkabilityScore, {}, 'C4 Walkability score', true);
Map.addLayer(C5_LastMileConnectivityScore, {}, 'C5 Last-mile score', false);
Map.addLayer(C6_PriorityInterventions, {}, 'C6 Priority interventions', false);
Map.addLayer(C1_UHI_HeatRisk, {}, 'C1 Heat risk', false);
Map.addLayer(C2_GreenShadeIndex, {}, 'C2 NDVI green/shade', false);
Map.addLayer(C3_TrafficRiskProxy, {}, 'C3 Traffic risk proxy', false);

// -------------------------
// 10) Exports to Drive
// -------------------------
var folder = 'GEE_Exports';
var scale = 10;

function exportMap(img, desc, prefix) {
  Export.image.toDrive({
    image: img,
    description: desc,
    folder: folder,
    fileNamePrefix: prefix,
    region: studyArea,
    scale: scale,
    maxPixels: 1e13
  });
}

exportMap(C1_UHI_HeatRisk, 'Indore_C1_HeatRisk', 'indore_c1_heat_risk');
exportMap(C2_GreenShadeIndex, 'Indore_C2_GreenShade_NDVI', 'indore_c2_green_shade_ndvi');
exportMap(C3_TrafficRiskProxy, 'Indore_C3_TrafficRisk_Proxy', 'indore_c3_traffic_risk_proxy');
exportMap(C4_WalkabilityScore, 'Indore_C4_WalkabilityScore', 'indore_c4_walkability_score');
exportMap(C5_LastMileConnectivityScore, 'Indore_C5_LastMileScore', 'indore_c5_lastmile_score');
exportMap(C6_PriorityInterventions, 'Indore_C6_PriorityInterventions', 'indore_c6_priority_interventions');
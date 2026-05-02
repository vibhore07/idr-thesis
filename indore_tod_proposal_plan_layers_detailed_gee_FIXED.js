/****
DETAILED TOD PROPOSAL PLAN (proposal layers only) for Indore Metro Stations
NO cumulativeCost used.

EXPORTS (GeoTIFF):
P1_TOD_ZonePlan_400_800_1500
P2_NodeHierarchy_and_Hubs
P3_Density_and_Height_Bands
P4_LandUseIntent_Framework
P5_NMT_and_CompleteStreet_Spine
P6_Parking_Management_Zones
P7_GreenBlue_PublicRealm_Framework
P8_InclusiveHousing_and_SocialInfra_Zones
P9_Implementation_Phasing_Map
****/

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

var r_core = 400, r_trans = 800, r_infl = 1500;
var pad = 7500;
var studyArea = stations.geometry().buffer(pad).bounds();
Map.centerObject(studyArea, 13);

// Base imagery
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(studyArea)
  .filterDate('2025-10-01', '2026-03-31')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .median()
  .clip(studyArea);
var s2Vis = {bands: ['B4','B3','B2'], min: 0, max: 3000, gamma: 1.15};

// WorldCover
var wc = ee.ImageCollection('ESA/WorldCover/v200').first().clip(studyArea);
var built = wc.eq(50).selfMask();
var green = wc.eq(10).or(wc.eq(20)).or(wc.eq(30)).or(wc.eq(40)).selfMask();
var water = wc.eq(80).selfMask();
var greenBlue = green.unmask(0).add(water.unmask(0)).gt(0).selfMask();

var builtDensity = built.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.circle({radius: 500, units: 'meters', normalize: true})
}).rename('built_density');

var distToGreen = greenBlue.fastDistanceTransform(30).sqrt().multiply(10).rename('dist_to_green_m');

// Heat (Landsat 9 LST)
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

var lstC = l9.select('ST_B10').multiply(0.00341802).add(149.0).subtract(273.15).rename('LST_C');
var lstSmooth = lstC.focal_mean({radius: 300, units: 'meters'}).rename('LST_C_smooth');

// Normalize
var bd = builtDensity.unitScale(0, 0.6).clamp(0, 1);
var gd = distToGreen.unitScale(0, 600).clamp(0, 1);
var ht = lstSmooth.unitScale(25, 45).clamp(0, 1);

// TOD zones
var coreZones = stations.map(function(f){ return f.buffer(r_core).copyProperties(f); });
var transZones = stations.map(function(f){ return f.buffer(r_trans).copyProperties(f); });
var inflZones = stations.map(function(f){ return f.buffer(r_infl).copyProperties(f); });

var g_core = coreZones.geometry().dissolve();
var g_trans = transZones.geometry().dissolve();
var g_infl = inflZones.geometry().dissolve();

var ring_core = g_core;
var ring_trans = g_trans.difference(g_core, 1);
var ring_infl = g_infl.difference(g_trans, 1);

var mask_all = ee.Image().byte().paint(g_infl, 1).selfMask();
var mask_core = ee.Image().byte().paint(ring_core, 1).selfMask();
var mask_trans = ee.Image().byte().paint(ring_trans, 1).selfMask();
var mask_infl = ee.Image().byte().paint(ring_infl, 1).selfMask();

// Station proximity (distance transform)
var stationSeeds = ee.Image().byte().paint(stations, 1).selfMask();
var distToStation = stationSeeds.fastDistanceTransform(30).sqrt().multiply(10).rename('dist_to_station_m');
var proxScore = ee.Image(1).subtract(distToStation.unitScale(0, r_infl).clamp(0, 1)).rename('prox').updateMask(mask_all);

// Hubs
var primaryStations = stations.filter(ee.Filter.eq('tier','primary'));
var secondaryStations = stations.filter(ee.Filter.eq('tier','secondary'));
var hubGeomPrimary = primaryStations.map(function(f){ return f.buffer(150); }).geometry().dissolve();
var hubGeomSecondary = secondaryStations.map(function(f){ return f.buffer(120); }).geometry().dissolve();

// Density bands
var densityBand = ee.Image(0)
  .where(mask_infl.mask(), 1)
  .where(mask_trans.mask(), 2)
  .where(mask_core.mask(), 3)
  .updateMask(mask_all);

var opportunityPocket = bd.lte(0.18).and(proxScore.gte(0.75)).and(mask_core.mask().or(mask_trans.mask())).selfMask();
var densityBand2 = densityBand.where(opportunityPocket, 3).rename('density_band2');

// Land use intent (conceptual)
var landIntent = ee.Image(0)
  .where(mask_infl.mask(), 4)
  .where(mask_trans.mask(), 3)
  .where(mask_core.mask(), 2);

var primaryHubMask = ee.Image().byte().paint(hubGeomPrimary, 1).selfMask();
landIntent = landIntent
  .where(bd.gte(0.35).and(primaryHubMask.mask()), 1)
  .where(greenBlue.unmask(0).gt(0), 5)
  .updateMask(mask_all)
  .rename('land_intent');

// Spine line connecting stations
var stationList = stations.toList(stations.size());
var coords = ee.List.sequence(0, stations.size().subtract(1)).map(function(i){
  return ee.Feature(stationList.get(i)).geometry().coordinates();
});
var spineLine = ee.Geometry.LineString(coords);
var completeStreetSpine = spineLine.buffer(40);
var nmtPriority = spineLine.buffer(150);

// Parking zones
var parkingZone = ee.Image(0)
  .where(mask_infl.mask(), 1)
  .where(mask_trans.mask(), 2)
  .where(mask_core.mask(), 3)
  .updateMask(mask_all)
  .rename('parking_zone');

// Green-blue & public realm
var coolingPriority = ht.gte(0.70).and(gd.gte(0.60)).and(mask_core.mask().or(mask_trans.mask())).selfMask();
var connectorPriority = gd.gte(0.60)
  .and(ee.Image().byte().paint(nmtPriority, 1).selfMask().mask())
  .and(mask_all.mask())
  .selfMask();

// Inclusive + social infra
var inclusiveZone = bd.lte(0.20)
  .and(proxScore.gte(0.55))
  .and(mask_trans.mask())
  .and(coolingPriority.mask().not())
  .selfMask();

var socialCatch = stations.map(function(f){ return f.buffer(600); }).geometry().dissolve();

// Phasing
var implScore = proxScore.multiply(0.40)
  .add(bd.multiply(0.25))
  .add(gd.multiply(0.20))
  .add(ht.multiply(0.15))
  .updateMask(mask_all);

var phase = ee.Image(3)
  .where(implScore.gte(0.66), 1)
  .where(implScore.gte(0.45).and(implScore.lt(0.66)), 2)
  .updateMask(mask_all)
  .rename('phase');

// Styling
function stationsStyled(fc) {
  var primary = fc.filter(ee.Filter.eq('tier','primary'))
    .style({color: 'd50000', pointSize: 10, width: 2, fillColor: 'd50000'});
  var secondary = fc.filter(ee.Filter.eq('tier','secondary'))
    .style({color: '2962ff', pointSize: 7, width: 2, fillColor: '2962ff'});
  return ee.ImageCollection([secondary, primary]).mosaic();
}
var stationsImg = stationsStyled(stations);
var boundaryImg = ee.Image().paint(studyArea, 1, 2).visualize({palette: ['000000'], opacity: 1});

var coreFill = ee.Image().paint(ring_core, 1).visualize({palette: ['ff6d00'], opacity: 0.22});
var transFill = ee.Image().paint(ring_trans, 1).visualize({palette: ['ffab40'], opacity: 0.15});
var inflFill = ee.Image().paint(ring_infl, 1).visualize({palette: ['ffe0b2'], opacity: 0.10});

var coreOutline = ee.Image().paint(ring_core, 1, 2).visualize({palette: ['ff6d00'], opacity: 1});
var transOutline = ee.Image().paint(ring_trans, 1, 2).visualize({palette: ['ffab40'], opacity: 1});
var inflOutline = ee.Image().paint(ring_infl, 1, 2).visualize({palette: ['ffe0b2'], opacity: 1});

var hubPrimImg = ee.Image().paint(hubGeomPrimary, 1).visualize({palette: ['b71c1c'], opacity: 0.35});
var hubSecImg  = ee.Image().paint(hubGeomSecondary, 1).visualize({palette: ['0d47a1'], opacity: 0.30});

var densityVis = densityBand2.visualize({min: 1, max: 3, palette: ['fff3e0','ffb74d','e65100'], opacity: 0.70});

var landVis = landIntent.visualize({
  min: 1, max: 5,
  palette: ['6a1b9a','d84315','fdd835','90a4ae','1b5e20'],
  opacity: 0.70
});

var completeStreetImg = ee.Image().paint(completeStreetSpine, 1).visualize({palette: ['263238'], opacity: 0.55});
var nmtPriorityImg = ee.Image().paint(nmtPriority, 1).visualize({palette: ['2e7d32'], opacity: 0.20});

var parkingVis = parkingZone.visualize({min: 1, max: 3, palette: ['cfd8dc','ffcc80','ff6d00'], opacity: 0.75});

var greenVis = green.visualize({palette: ['1b5e20'], opacity: 0.60});
var waterVis = water.visualize({palette: ['1565c0'], opacity: 0.75});
var coolingVis = coolingPriority.visualize({palette: ['ff1744'], opacity: 0.30});
var connectorVis = connectorPriority.visualize({palette: ['00c853'], opacity: 0.30});

var inclusiveVis = inclusiveZone.visualize({palette: ['7c4dff'], opacity: 0.35});
var socialCatchImg = ee.Image().paint(socialCatch, 1).visualize({palette: ['00b0ff'], opacity: 0.12});

var phaseVis = phase.visualize({min: 1, max: 3, palette: ['d32f2f','fbc02d','388e3c'], opacity: 0.75});

// Map products
var P1_TOD_ZonePlan_400_800_1500 = s2.visualize(s2Vis)
  .blend(inflFill).blend(transFill).blend(coreFill)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P2_NodeHierarchy_and_Hubs = s2.visualize(s2Vis)
  .blend(hubPrimImg).blend(hubSecImg)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P3_Density_and_Height_Bands = s2.visualize(s2Vis)
  .blend(densityVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P4_LandUseIntent_Framework = s2.visualize(s2Vis)
  .blend(landVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P5_NMT_and_CompleteStreet_Spine = s2.visualize(s2Vis)
  .blend(nmtPriorityImg).blend(completeStreetImg)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P6_Parking_Management_Zones = s2.visualize(s2Vis)
  .blend(parkingVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P7_GreenBlue_PublicRealm_Framework = s2.visualize(s2Vis)
  .blend(greenVis).blend(waterVis)
  .blend(connectorVis).blend(coolingVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P8_InclusiveHousing_and_SocialInfra_Zones = s2.visualize(s2Vis)
  .blend(socialCatchImg).blend(inclusiveVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

var P9_Implementation_Phasing_Map = s2.visualize(s2Vis)
  .blend(phaseVis)
  .blend(inflOutline).blend(transOutline).blend(coreOutline)
  .blend(stationsImg).blend(boundaryImg).clip(studyArea);

// Preview
Map.addLayer(P1_TOD_ZonePlan_400_800_1500, {}, 'EXPORT P1: TOD Zone Plan', true);
Map.addLayer(P2_NodeHierarchy_and_Hubs, {}, 'EXPORT P2: Nodes & Hubs', false);
Map.addLayer(P3_Density_and_Height_Bands, {}, 'EXPORT P3: Density/Height', false);
Map.addLayer(P4_LandUseIntent_Framework, {}, 'EXPORT P4: Land Use Intent', false);
Map.addLayer(P5_NMT_and_CompleteStreet_Spine, {}, 'EXPORT P5: NMT + Complete Street', false);
Map.addLayer(P6_Parking_Management_Zones, {}, 'EXPORT P6: Parking Zones', false);
Map.addLayer(P7_GreenBlue_PublicRealm_Framework, {}, 'EXPORT P7: Green-Blue & Public Realm', false);
Map.addLayer(P8_InclusiveHousing_and_SocialInfra_Zones, {}, 'EXPORT P8: Inclusive + Social', false);
Map.addLayer(P9_Implementation_Phasing_Map, {}, 'EXPORT P9: Phasing', false);

// Exports
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

exportMap(P1_TOD_ZonePlan_400_800_1500, 'Indore_TOD_P1_ZonePlan', 'indore_tod_p1_zone_plan');
exportMap(P2_NodeHierarchy_and_Hubs, 'Indore_TOD_P2_NodeHierarchy_Hubs', 'indore_tod_p2_nodes_hubs');
exportMap(P3_Density_and_Height_Bands, 'Indore_TOD_P3_Density_Height', 'indore_tod_p3_density_height');
exportMap(P4_LandUseIntent_Framework, 'Indore_TOD_P4_LandUseIntent', 'indore_tod_p4_landuse_intent');
exportMap(P5_NMT_and_CompleteStreet_Spine, 'Indore_TOD_P5_NMT_CompleteStreet', 'indore_tod_p5_nmt_complete_street');
exportMap(P6_Parking_Management_Zones, 'Indore_TOD_P6_ParkingZones', 'indore_tod_p6_parking_zones');
exportMap(P7_GreenBlue_PublicRealm_Framework, 'Indore_TOD_P7_GreenBlue_PublicRealm', 'indore_tod_p7_greenblue_publicrealm');
exportMap(P8_InclusiveHousing_and_SocialInfra_Zones, 'Indore_TOD_P8_Inclusive_SocialInfra', 'indore_tod_p8_inclusive_socialinfra');
exportMap(P9_Implementation_Phasing_Map, 'Indore_TOD_P9_Phasing', 'indore_tod_p9_phasing');
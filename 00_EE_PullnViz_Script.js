// This script is also available at 
// https://code.earthengine.google.com/f367c6b62933ef19a9c4036ea11c81d4


// load the required functions
var jrc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater")
var ls_fun = require('users/eeProject/water_surface_interpolation:functions_landsat.js');
var occ = jrc.select('occurrence');
var calcAreaImage = function(maskImage) {
  return(ee.Image.pixelArea().divide(1000000).mask(maskImage).rename(maskImage.bandNames()));
};
var createParagraphPanel = function(paragraph, labelStyle) {
  var labelList = [];
  for (var i = 0; i < paragraph.length; i++) {
    labelList.push(ui.Label(paragraph[i], labelStyle));
  } 
  return(ui.Panel(labelList));
}; 

// set parameters
var MAXDISTANCE = 1000; // used in cumulative_cost function to isolate lake area from water area
var waterClassification = 'Jones2019_3'; // one of ['Jones2019_3', 'Jones2019_2', 'Zou2018']


// load lake polygons
// load lake polygons
var lake1 = ee.FeatureCollection('users/sarinarinab/CitSci_Washington_Lakes').select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);
var lake2 = ee.FeatureCollection('users/sarinarinab/Horsepen');
var lake3 = ee.FeatureCollection('users/sarinarinab/NC_lakes');
var lake4 = ee.FeatureCollection('users/sarinarinab/CitSci_Washington_Lakes_newsummer2019').select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);
var lake5 = ee.FeatureCollection('users/sarinarinab/CitSci_Illinois_Lakes').select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);
var lake6 = ee.FeatureCollection('users/sarinarinab/Deep_Quarry_wName').select(['GNIS_Name1', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);
var lake7 = ee.FeatureCollection('users/sarinarinab/Herrick_Lake').select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);
var lake8 = ee.FeatureCollection('users/sarinarinab/CitSci_HarrierLakeIl').select(['name', 'areakm2'], ['GNIS_NAME', 'AREASQKM']);
//var lake9 = ee.FeatureCollection('users/sarinarinab/CitSci_FranceSHPfiles2').select(['name', 'areakm2'], ['GNIS_NAME', 'AREASQKM']);
//wisconsin lakes - can add back in
var lake10 = ee.FeatureCollection('users/sarinarinab/CitSci_Wisconsin_Lakes_20152019').select(['GNIS_Name', 'AreaSqKm'], ['GNIS_NAME', 'AREASQKM']);


var lakes = ee.FeatureCollection([lake1,lake2,lake3,lake4,lake5,lake6,lake7,lake8,lake10]).flatten()
var lakesNames = lakes.aggregate_array('GNIS_NAME');

// load landsat images
var dat = ls_fun.merge_collections_std_bandnames_collection1tier1_sr()
.filterDate('2015-01-01', '2030-01-01');
.filterMetadata('CLOUD_COVER', 'less_than', 30);   


///////////////////Where I've been working - Start///////////////////////////////////////////////////////////////////////////
//load sentinel 2 images
var s2 = ee.ImageCollection("COPERNICUS/S2_SR")
.filterDate('2015-01-01', '2030-01-01')
.filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', 30);   

//transform s2 images

var img = ee.Image(ee.ImageCollection(s2).first()).aside(print);  //something wrong w filterBounds - I think its cause S2 is an image collection - converted this line to be an image collection

var linearTransform = function(b, k, img, band) {
  return(img.select(band).multiply(k).add(b));
};

function cloudMaskS2(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(1).or(
             qa.bitwiseAnd(cirrusBitMask).eq(1));
             
  return(mask);
}

var transform = function(img) {
  
  var cloudmask = cloudMaskS2(img);
  // add snow/ice, and cloud shadow as value 1;
  cloudmask = cloudmask.or(img.select(['SCL']).eq(3).or(img.select(['SCL']).eq(11))).rename('obsMask');
  
  // Map.addLayer(img, {bands: ['B4', 'B3', 'B2'], min: 0, max: 10000, gamma: 1.5}, 'rgb');
  // Map.addLayer(qualityImg, {min: 1, max: 11, palette: ['ff0004', '868686', '774b0a', '10d22c', 'ffff52', '0000ff', '818181', 'c0c0c0', 'f1f1f1', 'bac5eb', '52fff9']}, 'classifiedImage');
  
  img = img.select(['B4', 'B3', 'B2', 'B8', 'B11', 'B12', 'QA60'], ['Red', 'Green', 'Blue', 'Nir', 'Swir1', 'Swir2', 'BQA']);  //select the bands and change the name of the cloud band - QA60 to BQA here ??
  
  var blue = linearTransform(0.0003, 0.9570, img, 'Blue');
  var green = linearTransform(0.0015, 1.0304, img, 'Green');
  var red = linearTransform(0.0041, 0.9533, img, 'Red');
  var nir = linearTransform(0.0139, 1.0157, img, 'Nir');
  var swir1 = linearTransform(0.0034, 0.9522, img, 'Swir1');
  var swir2 = linearTransform(0.0004, 0.9711, img, 'Swir2');
  
  var imgConverted = blue
  .addBands(green)
  .addBands(red)
  .addBands(nir)
  .addBands(swir1)
  .addBands(swir2)
  .addBands(cloudmask)
  .copyProperties(img)
  .copyProperties(img, ['system:index', 'system:time_start']);
  
  imgConverted = ee.Image(imgConverted);

  return(imgConverted);
};

//itierate over whole image collection - correct??
var s22 = s2.map(transform);

// print(transform(img));

// print(s22.first());
// print(dat.first());


// addCloudmask for landsat image
dat = dat.map(function(i) {
  var fmask = ls_fun.AddFmaskSR(i).select(['fmask']);
  return(i.select(['Blue', 'Green', 'Red', 'Swir1', 'Nir', 'Swir2']).addBands(fmask.gte(2).rename('obsMask')));
});
// 
// print(dat.first());


// merge landat and s2

var mergedCol = dat.merge(s22);

// print(mergedCol.first())

// print(ee.Image(mergedCol.filterMetadata('system:index', 'equals', '2_20181217T164709_20181217T164715_T16TCM').first().aside(print)).select(['Green']).projection().crs())


//////////////////////////////Where I've been working - End//////////////////////////////////////////////////////////////////



// get layers of the current map
var layers = Map.layers();

var lakeSelector = ui.Select(lakesNames.getInfo(), 'Select a lake');

lakeSelector.onChange(function(name) {
  
  var f = ee.Feature(lakes.filterMetadata('GNIS_NAME', 'equals', name).first());
  var aoi = f.geometry();
  
  Map.centerObject(aoi);
  
  var result = ee.ImageCollection(ls_fun.filterContained(mergedCol, aoi)) //changed from dat to mergedCol
  .map(function(img) {
    var i = ls_fun.CalculateWaterAddFlagsSR(img, waterClassification); // water classificaton method can be changed here
    var aoiBuffered = aoi.buffer(200);
    
    var waterMask = i.select(['waterMask']); // waterMask
    var obsMask = i.select(['obsMask']); // obstructionMask
    
    // reconstruct the water surface
    var recon = ls_fun.fillWater(waterMask, obsMask, aoiBuffered);
    
    // calculate lakeMask after reconstructing the water surface
    var reconstructedLake = ls_fun.ExtractChannel(recon.select(['Reconstructed']), aoi.buffer(-50), MAXDISTANCE)
    .rename(['Reconstructed']);
    var reconstructedLakeAll = ls_fun.ExtractChannel(recon.select(['ReconstructedAll']), aoi.buffer(-50), MAXDISTANCE)
    .rename(['ReconstructedAll']);
    var rawLake = waterMask.and(reconstructedLake).rename(['lakeMask']);
    recon = recon.select(['obsMask']).addBands(reconstructedLake).addBands(rawLake).addBands(reconstructedLakeAll);
    
    // construct area image, each pixel value = its area in km^2, each band masked with its corresponding lake mask
    var areaImage = calcAreaImage(recon.select(['lakeMask']))
    .addBands(calcAreaImage(recon.select(['obsMask'])))
    .addBands(calcAreaImage(recon.select(['Reconstructed'])))
    .addBands(calcAreaImage(recon.select(['ReconstructedAll'])));
    
    // add the pixel areas together to calculate total area
    var lakeArea = areaImage.reduceRegion({
      reducer: ee.Reducer.sum(), 
      geometry: aoiBuffered, 
      scale: 30});
      
    return(i.addBands(recon).setMulti(lakeArea).copyProperties(recon).set('refArea', f.get('AREASQKM')));
  });
  
  result
  .sort('timestamp')
  .filterMetadata('fillStatus', 'equals', 'success')
  
  var halfArea = aoi.area(10).divide(2).divide(1000000);
  
  var chart = ui.Chart.feature.byFeature({
    features: result
    // .filterMetadata('fillStatus', 'equals', 'success')
    // .filterMetadata('lakeMask', 'greater_than', halfArea)
    .sort('timestamp'), 
    xProperty: 'timestamp', 
    yProperties: ['refArea', 'lakeMask', 'Reconstructed', 'ReconstructedAll']})
    .setChartType('ScatterChart')
    .setOptions({
      title: 'Changes in lake area: ' + name + ' (size of lake polygon: ' + String(f.get('AREASQKM').getInfo()) + ' km^2)',
      vAxis: {title: 'Lake area (km^2)'},
      hAxis: {title: 'Date'},
      curveType: 'none',
      lineWidth: 1,
      pointSize: 4,
      explorer: {
        keepInBounds: true,
        actions: ['dragToZoom', 'rightClickToReset']
      },
      series: [
        {color: 'black', labelInLegend: 'Lake polygon area', pointSize: 0},
        {color: 'blue', labelInLegend: 'Lake area (exposed)'},
        {color: 'red', labelInLegend: 'Lake area recon 1'},
        {color: 'green', labelInLegend: 'Lake area recon 2'}
      ]
    });
    
  chart.style({
    stretch: 'vertical',
    minHeight: '800px'
  });
  
  chart.onClick(function(x, y, z) {
    var f = ee.Image(result
    .filterMetadata('timestamp', 'equals', x)
    .filterMetadata(z, 'equals', y)
    .first());
    // print(f);
    var imgId = f.get('image_id');
    var img = ee.Image(mergedCol.filterMetadata('system:index', 'equals', imgId).first());
    // print(img);
    Map.centerObject(aoi);
    layers.set(0, ui.Map.Layer(img, {bands: ['Red', 'Green', 'Blue'], min: 0, max: 2500}, 'rgb image', true));
    layers.set(1, ui.Map.Layer(img.normalizedDifference(['Green', 'Swir1']), {min: -1, max: 1, palette: ['#d8b365','#f5f5f5','#5ab4ac']}, 'mndwi image', false));
    layers.set(2, ui.Map.Layer(f.select(['lakeMask']).selfMask(), {palette: ['blue']}, 'exposed lake area', true));
    layers.set(3, ui.Map.Layer(f.select(['Reconstructed']).selfMask(), {palette: ['red']}, 'lake area (fill in clouded area only)', false, 0.8));
    layers.set(4, ui.Map.Layer(f.select(['ReconstructedAll']).selfMask(), {palette: ['green']}, 'lake area (fill in entire lake area)', true, 0.8));
    layers.set(5, ui.Map.Layer(f.select(['obsMask']).selfMask(), {palette: ['yellow']}, 'obstructions', false, 0.5));
    layers.set(6, ui.Map.Layer(occ.selfMask(), {palette: ['#9e0142','#d53e4f','#f46d43','#fdae61','#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd','#5e4fa2'], min: 1, max: 100}, 'jrc_occ', false, 0.7));
    layers.set(7, ui.Map.Layer(lakes, {color: 'grey'}, 'lake polygons', true, 0.5));
    
    //diagnose water classification
    var wm = ls_fun.CalculateWaterAddFlagsSR(img, waterClassification).select('waterMask');
    layers.set(8, ui.Map.Layer(wm.selfMask(), {palette: ['blue']}, 'waterClassification', false, 0.5));
    });
  
  chartPanel.widgets().set(0, chart);
  chartPanel.widgets().set(1, paragraphPanel);
  
  Export.table.toDrive({
    collection: result, 
    description: name + '_lake_area', 
    folder: '', 
    fileNamePrefix: name + '_lake_area', 
    fileFormat: 'csv'});
});

Map.setOptions('satellite');
Map.add(lakeSelector);
var chartPanel = ui.Panel([], ui.Panel.Layout.Flow('vertical'), {position: 'bottom-right', stretch: 'vertical', maxHeight: '800px'});
var paragraph = [
  'Lake polygon area: reference lake areas (km^2) from the AREASQKM property of the lake polygon file.',
  'Lake area (exposed) km^2: area of the lake surface that is not obstructed',
  'Lake area recon 1 km^2: lake area exposed + estimated cloud-covered lake area (estimated from water occurrence map)',
  'Lake area recon 2 km^2: lake area estimated directly from water occurrence map (Pekel et al, 2017).'
  ];
var txtStyleMain = {position: 'top-center', stretch: 'horizontal', textAlign: 'left', padding: '2px 10px 2px 10px', margin: '0px'};
var paragraphPanel = createParagraphPanel(paragraph, txtStyleMain);

ui.root.insert(1, chartPanel);
ui.root.setLayout(ui.Panel.Layout.flow('vertical', false));


#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Citizen Science Functions
Created on Mon Jul 23 14:59:36 2018
@author: simontopp
"""

#Define functions and set global variables

#Landsat Cloud Filters

###These are functions for unpacking the bit quality assessment band for Landsat 
def Unpack(bitBand, startingBit, bitWidth):
  #unpacking bit bands
  #see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
  return (ee.Image(bitBand)\
  .rightShift(startingBit)\
  .bitwiseAnd(ee.Number(2).pow(ee.Number(bitWidth)).subtract(ee.Number(1)).int()))
  
def UnpackAll(bitBand, bitInfo):
  unpackedImage = ee.Image.cat([Unpack(bitBand, bitInfo[key][0], bitInfo[key][1]).rename([key]) for key in bitInfo])
  return unpackedImage
  
def lsCloudFilter(image)  
  if collection == 'SR':
    mission = ee.String(lsSample.get('SATELLITE')).split('_').get(1)
    bitAffected = {
    'Cloud': [5, 1],
    'CloudShadow': [3, 1],
    'CirrusConfidence': [8,2], #only present for L8, but bits aren't used in 5-7 so will just come up empty
    'SnowIceConfidence': [9, 2],  #Be aware that super turbid water sometimes flagged as ice
    }
  
  elif collection == 'TOA':
    mission = ee.String(lsSample.get('SPACECRAFT_ID')).split('_').get(1)
    bitAffected = {
    'Cloud': [4, 1],
    'CloudShadow': [7, 2],
    'CirrusConfidence': [11,2]
    'SnowIceConfidence': [9, 2]
    }
    
  
  #Select qa band
  qa = image.select('qa')
    
  #Create cloud, shadow, and ice mask. 
  
  #Upack quality band to identify clouds, cloud shadows, and cirrus shadows.
  #For SR collections, clouds and cloud shadows will be either 0 or 1.  For
  #TOA collection, cloud shadow will be 1,2, or 3, associated with low, medium, or high
  #confidence respectively.  The following code only removes high confidence cloud shadows and cirrus
  #clouds, but this can be changed to accomadate specific research goals/areas.
  
  qaUnpack = UnpackAll(qa, bitAffected)
  
  if collection == 'SR':
    mask = qaUnpack.select('Cloud').eq(1)\
    .Or(qaUnpack.select('CloudShadow').eq(1))\
    .Or(qaUnpack.select('CirrusConfidence').eq(3))\
    .Or(qaUnpack.select('SnowIceConfidence').eq(3)))
  elif collection == 'TOA':
    mask = qaUnpack.select('Cloud').eq(1)\
    .Or(qaUnpack.select('CloudShadow').eq(3))\
    .Or(qaUnpack.select('CirrusConfidence').eq(3))\
    .Or(qaUnpack.select('SnowIceConfidence').eq(3))
  
  cScore = mask.reduceRegion(ee.Reducer.sum(),polygon.buffer(-60), 30).get('cloud')
  return image.set({"cScore": cScore})


# Clip image to Lake
def clipImage(image):
    return image.clip(buff)

'''
#Simple Landsat cloudscore function 
def cloudScore(image):
    scored = ee.Algorithms.Landsat.simpleCloudScore(image)
    cloudy = scored.select('cloud').gte(20)
    cScore = ee.Number(cloudy.reduceRegion(ee.Reducer.sum(),polygon.buffer(-60), 30).get('cloud'))
    return image.set({"cScore": cScore}).set({"cfStrict" : 0})

#For use with MNDWI

def fmaskScore(image):
    fmaskClouds = image.select('cfmask').eq(4)  #fmask band equal to 4:cloud
    cScore = fmaskClouds.reduceRegion(ee.Reducer.sum(), polygon, 30).get('cfmask')
    noWater = image.select('cfmask').neq(1)  #fmask band equal to 4:cloud
    cfStrict = noWater.reduceRegion(ee.Reducer.sum(), polygon.centroid().buffer(200), 30).get('cfmask')
    return image.set({"cScore" : cScore}).set({"cfStrict" : cfStrict})


# For use only with fmask 
def fmaskScore(image):
    fmaskClouds = image.select('cfmask').neq(4)  #fmask band equal to 4:cloud
    if lakename == 'Lake Mattamuskeet W' or 'Lake Mattamuskeet E' :
        cScore = fmaskClouds.reduceRegion(ee.Reducer.sum(), polygon.buffer(-650), 30).get('cfmask')
    #elif lakename == 'Lake Mattamuskeet E':
        #cSum = fmaskClouds.reduceRegion(ee.Reducer.sum(), polygon.buffer(-650), 30)
    else:
        cScore = fmaskClouds.reduceRegion(ee.Reducer.sum(), polygon.buffer(-100), 30).get('cfmask')
    return image.set({"cScore" : cScore})
''' 
##Reduce Resolution of S2 to match Landsat
# def reproject(image):
#     LSProjection = ee.Image(collectionL5.first()).projection()
#     S230m = image.reduceResolution(ee.Reducer.mean()).reproject(LSProjection)
#     return S230m.copyProperties(image).set({"system:time_start": image.get("system:time_start")})
#     # Force the next reprojection to aggregate instead of resampling
    
def reproject(image):
    LSProjection = ee.Image(collectionL5.first()).projection()
    S230m = image.reduceResolution(ee.Reducer.mean()).reproject(LSProjection)
    return S230m.copyProperties(image).set({"system:time_start": image.get("system:time_start")})

'''
#Sentinel 2 Cloud Score based on Quality Assessment band
def cloudScoreS2(image):
    qa = image.select('QA60').int16()
    #Clouds and cirrus are bit 10 and 11 respectivly
    cloudBitMask = 2**10
    cirrusBitMask = 2**11
    #Determines cloudy pixels (0 indicating clear pixels).
    cloudMask = qa.bitwiseAnd(cloudBitMask).gt(0).And(qa.bitwiseAnd(cirrusBitMask).gt(0))
    #scored = cloudMask.updateMask(cloudMask)
    cScore = ee.Number(cloudMask.reduceRegion(ee.Reducer.sum(),polygon).get('QA60')) 
    return image.set({"cScore": cScore})

#Simplified function based on QA band
def cloudS2Test(image):
  #Opaque and cirrus cloud masks cause bits 10 and 11 in QA60 to be set,
  #so values less than 1024 are cloud-free
    clouds = image.select('QA60').gte(1024)
    cDilate = clouds.focal_max(5)
    cScore = ee.Number(cDilate.reduceRegion(ee.Reducer.sum(), polygon, 10).get("QA60"))  
    return image.set({"cScore": cScore})

#Empirical Cloud Score Function for S2

def rescale(image, exp, thresholds):
    return image.expression(exp, {"image": image}).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])

def cloudScoreSII(image):
#Compute several indicators of cloudyness and take the minimum of them.
    score = ee.Image(1)
  #Clouds are reasonably bright in the blue and cirrus bands.
    score = score.min(rescale(image, "image.B2", [0.1, 0.5]))
    score = score.min(rescale(image, "image.B1", [0.1, 0.3]))
    score = score.min(rescale(image, "image.B1 + image.B10", [0.15, 0.2]))
  
  #Clouds are reasonably bright in all visible bands.
    score = score.min(rescale(image, "image.B4 + image.B3 + image.B2", [0.2, 0.8]))
  #Clouds are moist
    ndmi = image.normalizedDifference(["B8","B11"])
    score=score.min(rescale(ndmi, 'image', [-0.1, 0.1]))
  
  #However, clouds are not snow.
    ndsi = image.normalizedDifference(["B3", "B11"])
    score= score.min(rescale(ndsi, 'image', [0.8, 0.6])) 
    score = score.multiply(100).byte()
    return image.addBands(score)

    
def cloudS2(image):
    t = image.select(['B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000)
    t = t.addBands(image.select('QA60'))
    out = ee.Image(t.copyProperties(image))
    #out = out.select(['QA60', 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12'],['QA60','cb', 'blue', 'green', 'red', 're1','re2','re3','nir', 'nir2', 'waterVapor', 'cirrus','swir1', 'swir2'])
    scored = cloudScoreSII(out)
    cloudy = scored.select('constant').gte(20)
    cScore = cloudy.reduceRegion(ee.Reducer.sum(), polygon, 10).get('constant')
    return image.set({'cScore': cScore})
'''

##### S2 Cloudscore Attempt 3 ####
def cloudScoreSII(image):
    cld = image.select('B2').multiply(image.select('B9')).divide(10000)
    cloud = cld.gt(125)
    scale = ee.Image.pixelArea().sqrt()
    buf = cloud.fastDistanceTransform(500, 'pixels').sqrt().multiply(scale)
    clouds = cloud.add(buf.lt(500))
    #clouds = clouds.updateMask(cloud.add(buf.lt(500)))
    cScore = cloud.reduceRegion(ee.Reducer.sum(), polygon.buffer(-90), 10).get('B2')
    return image.set({'cScore': cScore})


####All in MNDWI - Green and Mid Infrared#####

#######  Otsu Dynamic Thresholding Method    ####
def thresh(image):
    ndwi = image.normalizedDifference(['Green','Swir1'])
    histogram = ndwi.reduceRegion(ee.Reducer.histogram(50).combine('mean', None, True).combine('variance', None, True), polygon.buffer(1500), 100, None, None, True)
    histogram = histogram.get('nd_histogram')

    def otsu(histogram):
        counts = ee.Array(ee.Dictionary(histogram).get('histogram'))
        means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
        size = means.length().get([0])
        total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
        sums = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
        mean = sums.divide(total)
        indices = ee.List.sequence(1, size)
  
  ### Compute between sum of squares, where each mean partitions the data.
        def thresh2(i):
            aCounts = counts.slice(0, 0, i)
            aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
            aMeans = means.slice(0, 0, i);
            aMean = aMeans.multiply(aCounts).reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount)
            bCount = total.subtract(aCount)
            bMean = sums.subtract(aCount.multiply(aMean)).divide(bCount)
            return aCount.multiply(aMean.subtract(mean).pow(2)).add(bCount.multiply(bMean.subtract(mean).pow(2)))
        bss = indices.map(thresh2)
        return means.sort(bss).get([-1])

        #### Return the mean value corresponding to the maximum BSS.
    threshold = otsu(histogram)
    return image.addBands(ndwi).rename('ndwi').set({'threshold': threshold})

        #### Return the mean value corresponding to the maximum BSS

'''
def Fmndwi(image):
    if x < 2:  ##
        ndwi = image.normalizedDifference(['B3','B6']) #L8
    else:
        ndwi = image.normalizedDifference(['B2','B5']) #L5/L7
    #ndwi = image.normalizedDifference(['B3','B8']) #S2
    return image.addBands(ndwi)
    
def otsu(image):
    ndwi = ee.Image(image).select('nd')
    histogram = ndwi.reduceRegion(ee.Reducer.histogram(50).combine('mean', None, True).combine('variance', None, True),polygon.buffer(1500),100, None, None, True)
    histogram = histogram.get('nd_histogram')
    def thresh(histogram):
        counts = ee.Array(ee.Dictionary(histogram).get('histogram'))
        means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'))
        size = means.length().get([0])
        total = counts.reduce(ee.Reducer.sum(), [0]).get([0])
        sums = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])
        mean = sums.divide(total)
        indices = ee.List.sequence(1, size)
  
  ### Compute between sum of squares, where each mean partitions the data.
        def thresh2(i):
            aCounts = counts.slice(0, 0, i)
            aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
            aMeans = means.slice(0, 0, i);
            aMean = aMeans.multiply(aCounts).reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount)
            bCount = total.subtract(aCount)
            bMean = sums.subtract(aCount.multiply(aMean)).divide(bCount)
            return aCount.multiply(aMean.subtract(mean).pow(2)).add(bCount.multiply(bMean.subtract(mean).pow(2)))
        bss = indices.map(thresh2)
        return means.sort(bss).get([-1])
    threshold = thresh(histogram)
        #### Return the mean value corresponding to the maximum BSS.
    return image.set({'threshold': threshold})



// Compute the histogram of the NDVI band.
var histogram = ndwi.reduceRegion({
  reducer: ee.Reducer.histogram(50)
      .combine('mean', null, true)
      .combine('variance', null, true), 
  // geometry: polygon, 
  scale: 100,
  bestEffort: true
});

print('histogram', histogram);

// main Otsu Function 
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  print(ui.Chart.array.values(ee.Array(bss), 0, means));
  
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

var threshold = otsu(histogram.get('nd_histogram')); 
print('threshold', threshold);

var classA = ndwi.gt(threshold);
Map.addLayer(classA.mask(classA), {palette: 'blue'}, 'class A');
'''


def waterExtentLS(image):
    threshold = ee.Number(image.get('threshold'))
    water = image.select("ndwi").gt(threshold)
    if lakename == 'Lake Mattamuskeet W' or 'Lake Mattamuskeet E':
        road = ee.FeatureCollection("TIGER/2016/Roads").filter(ee.Filter.eq('fullname', 'Frying Pan Lndg')).geometry().buffer(60)
        cost = water.eq(0).paint(road, 1)
    else:
        cost = water.eq(0)
    sources = ee.Image().toByte().paint(polygon.centroid().buffer(50), 1)
    sources = sources.updateMask(sources)
    cumulativeCost = cost.cumulativeCost(sources, 8000)
    area = cumulativeCost.lt(1).updateMask(area).multiply(ee.Image.pixelArea())
    #area = area.updateMask(area)
    #pArea = area.multiply(ee.Image.pixelArea())
    areaSum = area.reduceRegion(ee.Reducer.sum(),polygon.buffer(2000), 30) #need to use same resolution as input (30 landsat)
    return image.set({"Area": areaSum.get("cumulative_cost")})


def waterExtentS2(image):
    threshold = ee.Number(image.get('threshold'))
    #LSProjection = ee.Image(collectionL5.first()).projection()
    water = image.select("nd").gt(threshold)
    #water = ndwi.gt(threshold)
    if lakename == 'Lake Mattamuskeet W' or 'Lake Mattamuskeet E':
        road = ee.FeatureCollection("TIGER/2016/Roads").filter(ee.Filter.eq('fullname', 'Frying Pan Lndg')).geometry().buffer(60)
        cost = water.eq(0).paint(road, 1)
    else:
        cost = water.eq(0)
    sources = ee.Image().toByte().paint(polygon.centroid().buffer(50), 1)
    sources = sources.updateMask(sources)
    cumulativeCost = cost.cumulativeCost(sources, 8000)
    area = cumulativeCost.select('cumulative_cost').eq(0)
    area = area.updateMask(area)
    pArea = area.multiply(ee.Image.pixelArea())
    areaSum = pArea.reduceRegion(ee.Reducer.sum(),polygon.buffer(2000), 20)
    return image.set({"Area": areaSum.get('cumulative_cost')})

#after ndwi with Sentinel 2, call projection().nominalscale() to find out what the resolution is after normalized diff.


 
def LS_Selection(image):
    CleanLS = ee.Feature(lakes.get(i))
    return CleanLS.set({'Date':image.get('DATE_ACQUIRED')})\
    .set({"Area":ee.Number(image.get("Area"))})\
    .set({'cScore':ee.Number(image.get('cScore'))})\
    .set({'CloudCover':image.get('CLOUD_COVER')})\
    .set({'Mission': missions[x]})\
    .set({'cfStrict':ee.Number(image.get('cfStrict'))})\
    .set({'ID': image.get('system:index')})

    
def LS_SelectionSR(image):
    CleanLS = ee.Feature(lakes.get(i))
    return CleanLS.set({'Date':ee.Date(image.get('system:time_start'))})\
    .set({"Area":ee.Number(image.get("Area"))})\
    .set({'cScore':ee.Number(image.get('cScore'))})\
    .set({'CloudCover':image.get('CLOUD_COVER')})\
    .set({'Mission': missions[x]})\
    .set({'cfStrict':ee.Number(image.get('cfStrict'))})\
    .set({'ID': image.get('system:index')})
        
def S2_Selection(image):
    CleanS2 = ee.Feature(lakes.get(i))
    return CleanS2.set({'Date':ee.Date(image.get('system:time_start'))})\
    .set({"Area":ee.Number(image.get("Area"))})\
    .set({'cScore':ee.Number(image.get('cScore'))})\
    .set({'CloudCover':image.get('CLOUD_COVER')})\
    .set({'Mission': 'S2'})\
    .set({'cfStrict': 0})\
    .set({'ID': image.get('system:index')})
        

def getClosestTile(feature):  #Extracts path and row to filter feature collection and limit edge effects.
    def getDistance(f):
        return f.set({'distance': f.geometry().centroid().distance(pt, 10)})

    pt = feature.centroid().geometry()
    # Load the secondary collection: WRS-2 polygons.
    wrs = ee.FeatureCollection('ft:1_RZgjlcqixp-L9hyS6NYGqLaKOlnhSC35AB5M5Ll')

    # Define a spatial filter as geometries that intersect.
    spatialFilter = (ee.Filter.intersects(
    leftField = '.geo',
    rightField = '.geo',
    maxError = 10
    ))

    # Define a save all join.
    saveAllJoin = ee.Join.saveAll(matchesKey = 'tiles')

    poi = ee.FeatureCollection(ee.Feature(pt, {}))

    # Apply the join.
    intersectJoined = saveAllJoin.apply(poi, wrs, spatialFilter)
    tiles = ee.FeatureCollection(ee.List(ee.Feature(intersectJoined.first()).get('tiles')))

    tilesDist = tiles.map(getDistance).sort('distance', True)

    closestTile = ee.Feature(tilesDist.first())
    path = closestTile.get('PATH')
    row = closestTile.get('ROW')
    return feature.set({'path':path}).set({'row':row})

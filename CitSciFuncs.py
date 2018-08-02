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
  
def lsCloudFilter(image):  
  
  #Select qa band
  qa = image.select('qa')
    
  #Create cloud, shadow, and ice mask. 
  
  #Upack quality band to identify clouds, cloud shadows, and cirrus shadows.
  #For SR collections, clouds and cloud shadows will be either 0 or 1.  For
  #TOA collection, cloud shadow will be 1,2, or 3, associated with low, medium, or high
  #confidence respectively.  The following code only removes high confidence cloud shadows and cirrus
  #clouds, but this can be changed to accomadate specific research goals/areas.
  
  qaUnpack = UnpackAll(qa, bitAffected)
  
  if series == 'SR':
    mask = qaUnpack.select('Cloud').eq(1)\
    .Or(qaUnpack.select('CloudShadow').eq(1))\
    .Or(qaUnpack.select('CirrusConfidence').eq(3))\
    .Or(qaUnpack.select('SnowIceConfidence').eq(3))
  elif series == 'TOA':
    mask = qaUnpack.select('Cloud').eq(1)\
    .Or(qaUnpack.select('CloudShadow').eq(3))\
    .Or(qaUnpack.select('CirrusConfidence').eq(3))\
    .Or(qaUnpack.select('SnowIceConfidence').eq(3))
  elif series == 'S2':
    mask = qaUnpack.select('Cloud').eq(1)\
    .Or(qaUnpack.select('CirrusConfidence').eq(1))
  
  cScore = mask.reduceRegion(ee.Reducer.sum(),polygon.buffer(-60), 30).get('Cloud')
  return image.set({"cScore": cScore})

##Additional Cloud Score algorithm for Sentinel 2 Since QA band is mediocre.  Adapted from
#Chris Hewig and Ian Housman

def scale(image):
  t = image.select([ 'B1','B2','B3','B4','B5','B6','B7','B8','B8A', 'B9','B10', 'B11','B12']).divide(10000) #Rescale to 0-1
  t = t.addBands(image.select(['QA60']))
  out = ee.Image(t.copyProperties(image).copyProperties(image,['system:time_start']))
  return out.select(bnS2, bns)

def rescale(img, thresholds):
  """
  Linear stretch of image between two threshold values.
  """
  return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])
    
def sentinelCloudScore(img):
  """
  Computes spectral indices of cloudyness and take the minimum of them.
  
  Each spectral index is fairly lenient because the group minimum 
  is a somewhat stringent comparison policy. side note -> this seems like a job for machine learning :)
    
  originally written by Matt Hancher for Landsat imagery
  adapted to Sentinel by Chris Hewig and Ian Housman
  """
    
  # cloud until proven otherwise
  score = ee.Image(1)

  # clouds are reasonably bright
  score = score.min(rescale(img.select(['Blue']), [0.1, 0.5]))
  score = score.min(rescale(img.select(['Aerosol']), [0.1, 0.3]))
  score = score.min(rescale(img.select(['Aerosol']).add(img.select(['Cirrus'])), [0.15, 0.2]))
  score = score.min(rescale(img.select(['Red']).add(img.select(['Green'])).add(img.select('Blue')), [0.2, 0.8]))
  # clouds are moist
  ndmi = img.normalizedDifference(['Red4','Swir1'])
  score=score.min(rescale(ndmi, [-0.1, 0.1]))
  # clouds are not snow.
  ndsi = img.normalizedDifference(['Green', 'Swir1'])
  score=score.min(rescale(ndsi, [0.8, 0.6])).rename(['cloudScore'])
  
  score = score.multiply(100).byte()
  
  ##Set threshold, should be between 10-30 with 10 being more strict
  cScoreII = score.select(['cloudScore']).gt(15)\
  .reduceRegion(ee.Reducer.sum(),polygon.buffer(-60), 30).values().get(0)
  
  return img.addBands(score).set({'cScoreII':cScoreII})


def shadowMask(img,cloudMaskType):
    """
    Finds cloud shadows in images
    
    Originally by Gennadii Donchyts, adapted by Ian Housman
    """
    
    def potentialShadow(cloudHeight):
        """
        Finds potential shadow areas from array of cloud heights
        
        returns an image stack (i.e. list of images) 
        """
        cloudHeight = ee.Number(cloudHeight)
        
        # shadow vector length
        shadowVector = zenith.tan().multiply(cloudHeight)
        
        # x and y components of shadow vector length
        x = azimuth.cos().multiply(shadowVector).divide(nominalScale).round()
        y = azimuth.sin().multiply(shadowVector).divide(nominalScale).round()
        
        # affine translation of clouds
        cloudShift = cloudMask.changeProj(cloudMask.projection(), cloudMask.projection().translate(x, y)) # could incorporate shadow stretch?
        
        return cloudShift
  
    # select a cloud mask
    cloudMask = img.select(cloudMaskType)
    
    # make sure it is binary (i.e. apply threshold to cloud score)
    cloudScoreThreshold = 10
    cloudMask = cloudMask.gt(cloudScoreThreshold)

    # solar geometry (radians)
    azimuth = ee.Number(img.get('solar_azimuth')).multiply(math.pi).divide(180.0).add(ee.Number(0.5).multiply(math.pi))
    zenith  = ee.Number(0.5).multiply(math.pi ).subtract(ee.Number(img.get('solar_zenith')).multiply(math.pi).divide(180.0))

    # find potential shadow areas based on cloud and solar geometry
    nominalScale = cloudMask.projection().nominalScale()
    cloudHeights = ee.List.sequence(500,4000,500)        
    potentialShadowStack = cloudHeights.map(potentialShadow)
    potentialShadow = ee.ImageCollection.fromImages(potentialShadowStack).max()

    # shadows are not clouds
    potentialShadow = potentialShadow.And(cloudMask.Not())

    # (modified) dark pixel detection 
    darkPixels = toa.normalizedDifference(['Green', 'Swir2']).gt(0.3)

    # shadows are dark
    shadows = potentialShadow.And(darkPixels).rename(['shadows'])
    
    # might be scope for one last check here. Dark surfaces (e.g. water, basalt, etc.) cause shadow commission errors.
    # perhaps using a NDWI (e.g. green and nir)

    return img.addBands(shadows)


# Clip image to Lake
def clipImage(image):
  geoclip = ee.Geometry(image.geometry()).buffer(-5000)
  clip = image.clip(geoclip)
  return clip

#Function for filtering out partial images.  If returns true (1), 
#there are no empty pixels in image
def filterToNonEmpty(images, region, scale):
  def anyNonZero(image):
    nonZero = image.select(0).reduceRegion(ee.Reducer.anyNonZero(),region.centroid().buffer(60),scale).values().get(0)
    return image.set({'anyNonZero': nonZero})
  return images.map(anyNonZero).filter(ee.Filter.eq('anyNonZero', 1))

###anyNonZero test
def filterToNonEmpty1(images, region, scale):
  def anyNonZero(image):
    nonZero = image.select(0).reduceRegion(ee.Reducer.allNonZero(),region,scale).values().get(0)
    return image.set({'anyNonZero': nonZero})
  return images.map(anyNonZero).filter(ee.Filter.eq('anyNonZero', 1))

def lakeClip(image):
  return image.clip(buff)

def reproject(image):
    S230m = image.reduceResolution(ee.Reducer.mean()).reproject( 'EPSG:32628', None, 30)
    return S230m.copyProperties(image).set({"system:time_start": image.get("system:time_start")})

####All in MNDWI - Green and Mid Infrared#####

#######  Otsu Dynamic Thresholding Method    ####
def thresh(image):
    if series == 'S2':
      ndwi = image.normalizedDifference(['Green','Swir1']).rename('ndwi')
      ndwi = ndwi.reproject(ndwi.projection().scale(3,3))
    else:
      ndwi = image.normalizedDifference(['Green','Swir1']).rename('ndwi')
    histogram = ndwi.reduceRegion(ee.Reducer.histogram(50).combine('mean', None, True).combine('variance', None, True), polygon.buffer(1500), 100, None, None, True)
    histogram = histogram.get('ndwi_histogram')

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
    return image.addBands(ndwi).set({'threshold': threshold})

        #### Return the mean value corresponding to the maximum BSS

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
    area = cumulativeCost.lt(1)
    area = area.updateMask(area).multiply(ee.Image.pixelArea())
    #area = area.updateMask(area)
    #pArea = area.multiply(ee.Image.pixelArea())
    areaSum = area.reduceRegion(ee.Reducer.sum(),polygon.buffer(2000), 30) #need to use same resolution as input (30 landsat)
    return image.set({"Area": areaSum.get("cumulative_cost")})

def out(image):
    CleanLS = ee.Feature(lakes.get(i))
    return CleanLS.set({'Date':ee.Date(image.get('system:time_start'))})\
    .set({"Area":ee.Number(image.get("Area"))})\
    .set({'cScore':ee.Number(image.get('cScore'))})\
    .set({'CloudCover':image.get('CLOUD_COVER')})\
    .set({'ID': image.get('system:index')})

        
#Extracts path and row to filter feature collection and limit edge effects.
#Compliments of Xiao

##Function for limiting the max number of tasks sent to
#earth engine at one time to avoid time out errors
#Compliments of Xiao

def maximum_no_of_tasks(MaxNActive, waitingPeriod):
  ##maintain a maximum number of active tasks
  time.sleep(10)
  ## initialize submitting jobs
  ts = list(ee.batch.Task.list())

  NActive = 0
  for task in ts:
       if ('RUNNING' in str(task) or 'READY' in str(task)):
           NActive += 1
  ## wait if the number of current active tasks reach the maximum number
  ## defined in MaxNActive
  while (NActive >= MaxNActive):
      time.sleep(waitingPeriod) # if reach or over maximum no. of active tasks, wait for 2min and check again
      ts = list(ee.batch.Task.list())
      NActive = 0
      for task in ts:
        if ('RUNNING' in str(task) or 'READY' in str(task)):
          NActive += 1
  return()    

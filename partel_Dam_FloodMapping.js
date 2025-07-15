//definin our region of intrest
var roi = ee.Geometry.Polygon([
  [
    [36.13374119328197, -0.10503321478618759],
    [36.15004902409252, -0.10503321478618759],
    [36.15004902409252, -0.09340317395527557],
    [36.13374119328197, -0.09340317395527557],
    [36.13374119328197, -0.10503321478618759]  
  ]
]);

Map.centerObject(roi,12)
Map.addLayer(roi,{},"ROI")

var collection = ee.ImageCollection('COPERNICUS/S1_GRD')
        .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
        .filter(ee.Filter.eq('instrumentMode', 'IW'))
        .filter(ee.Filter.or(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'),ee.Filter.eq('orbitProperties_pass', 'ASCENDING')))
        
        
var before = collection.filter(ee.Filter.date('2018-04-01', '2018-05-07')).filterBounds(roi)    
var after = collection.filter(ee.Filter.date('2018-05-10', '2018-05-25')).filterBounds(roi)

print(before)
print(after)

var before_image = before.select('VH').mosaic().clip(roi);
var after_image = after.select('VH').mosaic().clip(roi);

var before_filtered = ee.Image(toDB(RefinedLee(toNatural(before_image))))
var after_filtered = ee.Image(toDB(RefinedLee(toNatural(after_image))))

var flood = before_filtered.gt(-20).and(after_filtered.lt(-20));
var flood_mask = flood.updateMask(flood.eq(1));

var water = before_filtered.lt(-20).and(after_filtered.lt(-20));
var water_mask = water.updateMask(water.eq(1));

Map.centerObject(roi);
//Map.addLayer(before_image,{min:-25,max:0},'before')
//Map.addLayer(after_image,{min:-25,max:0},'after')
Map.addLayer(before_filtered,{min:-25,max:0},'before_filtered')
Map.addLayer(after_filtered,{min:-25,max:0},'after_filtered')
// Map.addLayer(flood_mask,{palette:['Yellow']},'Flood_Inundation')
Map.addLayer(water_mask,{palette:['Blue']},'Water')

print('Total ROI Area (Ha)', roi.area(1).divide(10000))

var stats = flood_mask.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: roi,
  scale: 10,
  maxPixels: 1e13,
  tileScale: 16
})

print(stats);
var flood_area = ee.Number(stats.get('sum')).divide(10000).round();
print('Flooded Area (Ha)', flood_area)

// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
// Create legend title
var legendTitle = ui.Label({
  value: 'Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette =['Yellow', 'Blue'];
 
// name of the legend
var names = ['Flood Inundation','Water body'];
 
// Add color and and names
for (var i = 0; i <2; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);

// set position of panel
var title = ui.Panel({
  style: {
    position: 'top-center',
    padding: '8px 15px'
  }
});
 
// Create legend title
var mapTitle = ui.Label({
  value: 'PARTEL DAM BURST FLOOD MAPPING 2018',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
title.add(mapTitle);
Map.add(title);


// Speckle Filtering Functions

// Function to convert from d
function toNatural(img) {
  return ee.Image(10.0).pow(img.select(0).divide(10.0));
}

//Function to convert to dB
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0);
}

//Apllying a Refined Lee Speckle filter as coded in the SNAP 3.0 S1TBX:

//https://github.com/senbox-org/s1tbx/blob/master/s1tbx-op-sar-processing/src/main/java/org/esa/s1tbx/sar/gpf/filtering/SpeckleFilters/RefinedLee.java
//Adapted by Guido Lemoine

// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return(result.arrayFlatten([['sum']]));
}

// // Exporting the masked flood_mask
// Export.image.toDrive({
//   image: flood_mask,
//   description: 'flood_mask',
//   folder: 'Microwave_RS',
//   fileNamePrefix: 'flooded_area',
//   region: roi,
//   scale: 10,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });


              //create an RGB image using sentinel 2
var s2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
                  .select(['B4','B3','B2']);

var rgbVis = {
  min: 0,
  max: 3000,
  bands: ['B4', 'B3', 'B2'],
};

        
var s2_before = s2.filter(ee.Filter.date('2018-03-10', '2018-03-25')).filterBounds(roi).mosaic().clip(roi)  
var s2_after = s2.filter(ee.Filter.date('2018-06-10', '2018-06-25')).filterBounds(roi).mosaic().clip(roi)

// print(before.size())
// print(after.size())


Map.addLayer(s2_after,rgbVis,"s2_after")
Map.addLayer(s2_before,rgbVis,"s2_before")


// // exporting the before and after RGB sentinel 2 data
// Export.image.toDrive({
//   image: s2_before,
//   description: 's2_before_flooding',
//   folder: 'Microwave_RS',
//   fileNamePrefix: 's2_before',
//   region: roi,
//   scale: 10,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });

// Export.image.toDrive({
//   image: s2_after,
//   description: 's2_after_floods',
//   folder: 'Microwave_RS',
//   fileNamePrefix: 's2_after',
//   region: roi,
//   scale: 10,
//   crs: 'EPSG:4326',
//   maxPixels: 1e13
// });

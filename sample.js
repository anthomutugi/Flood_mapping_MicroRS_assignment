// Displaying and visualizing the image on the map interface
Map.centerObject(image,12)
Map.addLayer(image)

// var taitataveta=ee.filter(ee.Filter.eq('COUNTY_NAM',TAITATAVETA));
// Displaying and visualizing the aoi on the map interface
Map.addLayer(aoi)

//Clipping the image within my aoi
Map.addLayer(image.clip(aoi),imageVisParam,'image1')
// Function to fill in the scanline gaps landsat 7
function fillGaps(image1) {
    var filled1 = image1.focal_mean(1, 'square', 'pixels', 5);
    var filled2 = image1.focal_mean(2, 'square', 'pixels', 5);
    var filled3 = image1.focal_mean(3, 'square', 'pixels', 5);
    
  // Combine the different scales of filling to create a composite image that is filled
    var NYD = image1.unmask(filled1).unmask(filled2).unmask(filled3);
    
    return NYD;
    }
//   //Defining of the image and applying of gap filling function
    var NYD = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2") 
            .filterDate('2005-01-01', '2005-12-31')
            .filterBounds(aoi)
            .filter(ee.Filter.lt('CLOUD_COVER',5))
            .map(fillGaps)//Application of gap filing function
            .median()
            .clip(aoi);
             
//   //visualization parameters to be used in the gap filled image
  var vis_params = { 
    bands:('SR_B1','SR_B2','SR_B3','SR_B4','SR_B5'),
    max:3000,
    min:0,
    gamma: 1,
    
    };         
  Map.addLayer(NYD, vis_params ,'NYDgap_filled')
  
  var NYDgap_filled= NYD.select(['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5' ]);
             
    // Exporting the image
// Export.image.toDrive({
//   image: NYDgap_filled, 
//   description: 'Gap_filled_nyandarua', 
//   folder: 'DIAimages', 
//   fileNamePrefix: 'NYANDARUA COUNTY FROM 2005-01-01 TO 2005-12-31 ', 
//   region: aoi, 
//   scale: 30, 
//   crs: 'EPSG:4326', 
//   maxPixels: 1e10, 
//   });
    
    
    
    
    
    
    
    
    
    
    
    
    

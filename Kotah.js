var aster = ee.ImageCollection("ASTER/AST_L1T_003"),
    sentinel = ee.ImageCollection("COPERNICUS/S2_SR"),
    landsat = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT"),
    sent = ee.ImageCollection("COPERNICUS/S2_SR"),
    roi = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Point([72.13059794327728, 34.550499453749666]),
    region2 = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[72.09403406998626, 34.60618644705522],
          [72.09403406998626, 34.517691696251326],
          [72.27565180680267, 34.517691696251326],
          [72.27565180680267, 34.60618644705522]]], null, false);
 //var point = ee.Geometry.Point([ 72.408198, 34.296350]); 
//***************************** Mardan gabbro Mapping using Pricipal Component Analysis and Decorrelation Stretching **************************************


//=================================================================== Read Me =========================================================///
/* _______________________________________THIS FILE CONTAINS ALL THE CODE USED FOR Granite MARDAN PAPER ________________________________

1) Sentinel image was filtered and experimented with to obtain the most suitable image which we loaded through image ID. The filtering code is commented
2) All the scatter plots/charts, bands TCCs, FCCs, and Individual PCs outputs are commented to avoid overloading the system. You can uncomment the 'print(...)' and Map.addLayer(...) to get these outputs
3) The code is ordered as follows:
  a) Filtering and Loading Images, Selecting bands
  b) NDVI calculation and obtaining binary NDVI map 
  c) DS and PCA main function definition
  d) Bands TCCs and FCCs (Outputs are commented)
  e) Apply DS algorithm on all dataset to stretch bands data
  f) Apply PCA of stretched bands data of Landsat -8, ASTER, Sentinel - 2 (in this order) output results and plot PC1 vs PC2 and PC3 scatter plots for each dataset
  g) Apply PCA on Raw data for Landsat -8, ASTER and Sentinel - 2 in this order.
  h) Obtain weighted linear combination of three best components from these datasets
  i) Center the map to the specified coordinates (Last line of code)
4) To get invidual PCs copy and paste this code snippet below PCA application:
5) 

for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImage.bandNames().get(i).getInfo();
  //Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}

*/
//==================================================== Read Me ============================================================////


var sentImage = ee.ImageCollection('COPERNICUS/S2_SR')
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 5))
.filter(ee.Filter.lt('NOT_VEGETATED_PERCENTAGE', 60))
.filterDate("2019-12-02", "2020-12-30");
//print(sentImage)
var sentImage = sentImage.median().clip(region2);

//Select Sentinel - 2 Bands to use in the study
var sentbands = ['B2','B3','B4','B8','B11','B12'];

var NDVI = sentImage.expression(
  "(NIR-RED) / (NIR+RED)",
  {
    RED: sentImage.select("B4"),
    NIR: sentImage.select("B8"),
    BLUE: sentImage.select("B2")
  })
  
  Map.addLayer(NDVI, {min: 0, max: 1}, "NDVI");

var region = region2;
var DSregion = region2;

// Nice visualization parameters for a vegetation index.
var vis = {min: 0, max: 1, palette: [
  'FFFFFF', 'CE7E45', 'FCD163', '66A000', '207401',
  '056201', '004C00', '023B01', '012E01', '011301']};
  /*
Map.addLayer(NDVI, vis, 'NDVI');
//Map.addLayer(NDVI.gt(0.5), vis, 'NDVI Binarized');
*/
var landCollection = ee.Image(landsat
.filterDate("2015-01-01", "2020-01-01")
.filterBounds(region).filterMetadata("CLOUD_COVER", "less_than", 10).median().clip(region));
//print("Landsat 8 Scene", landCollection);

var landImage = landCollection.clip(region2);
var landbands = ['B5','B6','B7','B10','B11'];

var asterImage = ee.Image(aster
.filterDate("2001-01-01", "2007-10-01")
.sort('CLOUDCOVER')
.filterBounds(region)
.first());
//print("ASTER Scene",asterImage);
var asterImage = asterImage.clip(region2);
//Select ASTER Bands to use in the study
var asterbands = ['B01', 'B02', 'B3N', 'B10', 'B11', 'B12'];

//Ignore this, this is only used to avoid errors while concatinting eigen pairs
var eigenCollection = ee.Array([[0],[1],[2],[3],[4]]);
var tscale = 5;

// Decorrelation Stretching Main Function     
function decStr(bandsImage, location, scale){
  var bandNames = bandsImage.bandNames();
  // Naming the axis for intuition
  var dataAxis = 0;
  var bandsAxis = 1;
  // Calculate the mean for each band image
  var meansAll = bandsImage.reduceRegion(ee.Reducer.mean(), location, scale);
  // Generate an array (1D Matrix) of mean of each band
  var arrayOfmeans = ee.Image(meansAll.toArray());
  // Collapse the bands data such that each pixel is a matrix of pixel values of each band
  var pixelArrays = bandsImage.toArray();
  // Use the means array and the collapsed band data to center each pixel of each band by subtracting its corresponding mean from it
  var meanCent = pixelArrays.subtract(arrayOfmeans);
  // Calculate the Covariance matrix for the bands data
  var covar = meanCent.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale
  });
  
  // Get the covariance in array format which shows the band-band covarince of the data
  var covarArray = ee.Array(covar.get('array'));
  // Perform eigen decomposition of the covariance matrix to obtain eigen values and eigen vector pairs
  var eigenPairs = covarArray.eigen();
  var eigenValues = eigenPairs.slice(bandsAxis, 0, 1); // slice(axis, start, end, step)
  var eigenVectors = eigenPairs.slice(bandsAxis, 1);
  // Rotate by the eigenvectors, scale to a variance of 30, and rotate back.
  //Store a diagonal matrix in i
  var i = ee.Array.identity(bandNames.length()); // i will be used to isolate each band data and scale its variance e.g i = [1,0,0,0,0] = isolate first band from 5 bands
  // Calculate variance from the eigenvalues ---> variance = 1/sqrt(eigenvalues)
  // matrixToDiag = Computes a square diagonal matrix from a single column matrix for multiplication purposes
  var variance = eigenValues.sqrt().matrixToDiag();
  //Multiply diagonal matrix i by 30 and divide by vaiance to obtain scaling variance matrix
  var scaled = i.multiply(30).divide(variance); //Changed from 30 -> 50, It was observed that changing variance scale increases contrast. Best contrast obtained for 30
  // Calculate a rotation matrix ---> rotationMatrix =  Eigenvect.Transpose * ScaledVariance * Eigenvect
  var rotation = eigenVectors.transpose()
    .matrixMultiply(scaled)
    .matrixMultiply(eigenVectors);
  // Convert 1-D nomalized array image data to 2-D and transpose it so it can be multiplied with rotation matrix
  var transposed = meanCent.arrayRepeat(bandsAxis, 1).arrayTranspose();
  // Multiply the transposed data with the rotation matrix
  return transposed.matrixMultiply(ee.Image(rotation))
    .arrayProject([bandsAxis])   //This drop unecessary axis from the transposed data and only retains 2 axis
    .arrayFlatten([bandNames])  //Flatten collections of collections
    .add(127).byte(); // Conver pixel values to 127 means so it can be visualized between 0 - 255 range.
    
    // .byte is used to force element wise operation
}

// Principal Component Analysis Main Function
function PCA(meanCent, scale, location){
  // Flatten the band image data in from 2D to a 1D array
  var arrays = meanCent.toArray();
  //print('PCA applying on', meanCent);
  // Calculate the covariance matrix for the bands data of the region
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale,
    //tileScale: tscale,
    maxPixels: 1e9
  });
  // Get the band to band covariance of the region in 'array' format. Here .get('array') --> casts to an array
  var covarArray = ee.Array(covar.get('array'));
  // Perform an eigen analysis and slice apart the values and vectors.
  var eigenPairs = covarArray.eigen();
  // This is a P-length vector of Eigenvalues. Here P = number of PCs
  var eigenValues = eigenPairs.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigenPairs.slice(1, 1);
  //Print and store eigen pairs in eigenCollection variable and export to drive
  print('eigen Values', eigenValues);
  print('eigen Vector', eigenVectors);
    //Make feature collection out of eigenpairs so it can be exported to excel. From there we Convert it to a table using a python script
  eigenCollection = ee.Feature(null,{values:ee.Array.cat([eigenValues,eigenVectors],1)}); 
  print('Eigen Collection Length',eigenCollection);
    // Export the FeatureCollection to excel sheet in drive

  Export.table.toDrive({
  collection: ee.FeatureCollection([eigenCollection]),
  description: 'eigenAnalysis',
  fileFormat: 'CSV'
  });

  // Convert the 1D image array back to 2D matrix for multiplication
  var imageMat = arrays.toArray(1);
  // To obtain PC = EigenVectors * 2D Image Matrix
  var PCs = ee.Image(eigenVectors).matrixMultiply(imageMat);
  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
  // Turn the PCs into a P-band image, normalized by SD.
  return PCs
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
}

          //TCCs and FCCs
  //ASTER L1T FCC
var trueColor = {
  bands: ["B05", "B04", "B3N"],
  min: 0,
  max: 300
};
//Map.addLayer(asterImage, trueColor, "ASTER False-color");

  //Sentinel - 2 L2A TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 True Color");

  //Sentinel - 2 L2A FCC
var trueColor = {
  bands: ["B11", "B12", "B8"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 False Color");

  //Landsat - 8 raw TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 30000
};
//Map.addLayer(landImage, trueColor, "Landsat True Color");

  //Landsat - 8 raw FCC
var trueColor = {
  bands: ["B6", "B7", "B5"],
  min: 0,
  max: 30000
};
//Map.addLayer(landImage, trueColor, "Landsat False Color");





              //APPLYING DS ON ALL DATA
// Selecting bands to apply DS
var landBandsImage = landImage.select(landbands);
var asterBandsImage = asterImage.select(asterbands);
var sentBandsImage = sentImage.select(sentbands);

//Obtain DS Results for All Satelites using dcs function
var DSLand = decStr(landBandsImage, region2, 1000);
var DSaster = decStr(asterBandsImage, region2, 1000);
var DSsent = decStr(sentBandsImage, region2, 1000);

//FCC of 3 bands of DS results for all satelites
var selectBands = [0,1,2,3,4];
//Map.addLayer(DSLand.select(selectBands), {}, 'DCS Landsat Image');
var selectBands = [0,1,2,3,4,5]; 
//Map.addLayer(DSaster.select(selectBands), {}, 'DCS Aster Image');
var selectBands = [0,1,2,3,4,5]; 
//Map.addLayer(DSsent.select(selectBands), {}, 'DCS Sentinel Image');




            //PRINCIPAL COMPONENT ANALYSIS
      //Applying PCA on DS of Landsat 8
// Obtain the geometry of the region from the raw bands since stretched bands doesnt have geometry
var region = landImage.geometry();
var bands = [0,1,2,3,4];
var image =  DSLand.select(bands);
// Set some paramters for PCA application.
var scale = 30;
var bandNames = image.bandNames();
//Obtain means and stored it in a dictionary
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
// For faster covariance reduction we mean center the data
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
// This function returns the names of bands in list form
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
// Here wer used the centered imagery and the parameters defined to obtain PC of the region specified
var pcImageLDS = PCA(centered, scale, region);
 /*for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImageLDS.bandNames().get(i).getInfo();
  Map.addLayer(pcImageLDS.select([band]), {min: -2, max: 2}, band);
}*/
/*Export.image.toDrive({
  image: pcImage,
  description: 'Landsat8PCAofDS',
  folder: "GraniteExports",
  scale: 30,
  region: region
});*/
//print("Landsat8PCAofDS", pcImageLDS);
// Obtain FCC of the selected PCs based on Crosta Technique
//Map.addLayer(pcImageLDS, {bands: ['pc2', 'pc3', 'pc5'], min: -2, max: 2}, 'Landsat - 8 PCA of DS PC 2,3,5');



  //Applying PCA on DS of ASTER
var region = asterImage.geometry();
var bands = [0,1,2,3,4,5];
var image =  DSaster.select(bands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImageADS = PCA(centered, scale, region);
 /*for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImageADS.bandNames().get(i).getInfo();
  Map.addLayer(pcImageADS.select([band]), {min: -2, max: 2}, band);
}*/
/*Export.image.toDrive({
  image: pcImageADS,
  description: 'ASTERPCAofDS',
  folder: "GraniteExports",
  scale: 30,
  region: region
});*/ 
//Map.addLayer(pcImageADS, {bands: ['pc3', 'pc4', 'pc5'], min: -2, max: 2}, 'ASTER L1T - PCA of DS PC 3,4,5');



  //Applying PCA on DS of Sentinel - 2
var region = sentImage.geometry();
var bands = [0,1,2,3,4,5];
var image =  DSsent.select(bands);
var scale = 10;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImageSDS = PCA(centered, scale, region);
 /*for (var i = 0; i < bandNames.length().getInfo(); i++) {
   var band = pcImageSDS.bandNames().get(i).getInfo();
  Map.addLayer(pcImageSDS.select([band]), {min: -2, max: 2}, band);
}*/
/*Export.image.toDrive({
  image: pcImageDS,
  description: 'Sentinel2PCAofDS',
  folder: "GraniteExports",
  scale: 15,
  region: region
});*/

//Map.addLayer(pcImageSDS, {bands: ['pc2', 'pc3', 'pc6'], min: -2, max: 2}, 'Sentinel 2 L2C - PCA of DS (pc 2,3,6)'); //changed from PC1, PC2 and PC3




/*Export.image.toDrive({
  image: pcImage,
  description: 'Landsat8PCA',
  folder: "GraniteExports",
  scale: 30,
  region: region
});*/ 
//Map.addLayer(pcImage, {bands: ['pc2', 'pc3', 'pc1'], min: -2, max: 2}, 'Landsat 8 - PCA');




      //ASTER L1T PCA
  //PCA
var region = asterImage.geometry();
var image =  asterImage.select(asterbands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);

var pcImageA = PCA(centered, scale, region);
/*Export.image.toDrive({
  image: pcImage,
  description: 'ASTERL1tPCA',
  folder: "GraniteExports",
  scale: 30,
  region: region
});*/
Map.addLayer(pcImageA, {bands: ['pc2', 'pc3', 'pc6'], min: -2, max: 2}, 'ASTER L1T - PCA');

  // Plot each PC as a new layer
  /*
for (var i = 0; i < bandNames.length().getInfo(); i++) {
var band = pcImageA.bandNames().get(i).getInfo();
  Map.addLayer(pcImageA.select([band]), {min: -2, max: 2}, band);
}*/

 //Landsat PCA
var region = landImage.geometry();
var image =  landImage.select(landbands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);


var pcImageL = PCA(centered, scale, region);
/*
for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImageL.bandNames().get(i).getInfo();
  Map.addLayer(pcImageL.select([band]), {min: -2, max: 2}, band);
}*/
/*Export.image.toDrive({
  image: pcImage,
  description: 'Sentinel2PCA',
  folder: "GraniteExports",
  scale: 15,
  region: region
});*/
// Plot each PC as a new layer
Map.addLayer(pcImageL, {bands: ['pc1', 'pc2', 'pc4'], min: -2, max: 2}, 'Landsat 8 - PCA  used in paper');


      
  //Sentinel PCA
var region = sentImage.geometry();
var image =  sentImage.select(sentbands);
var scale = 20;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);


var pcImageS = PCA(centered, scale, region);
/*
 for (var i = 0; i < bandNames.length().getInfo(); i++) {
 var band = pcImageS.bandNames().get(i).getInfo();
  Map.addLayer(pcImageS.select([band]), {min: -2, max: 2}, band);
}*/
/*Export.image.toDrive({
  image: pcImage,
  description: 'Sentinel2PCA',
  folder: "GraniteExports",
  scale: 15,
  region: region
});*/
// Plot each PC as a new layerl
Map.addLayer(pcImageS, {bands: ['pc2', 'pc3', 'pc5'], min: -2, max: 2}, 'Sentinel 2 - PCA  124');
/*Map.addLayer(pcImageS, {bands: ['pc1', 'pc2', 'pc3'], min: -2, max: 2}, 'Sentinel 2 - PCA  123');
Map.addLayer(pcImageS, {bands: ['pc1', 'pc2', 'pc5'], min: -2, max: 2}, 'Sentinel 2 - PCA  125');
Map.addLayer(pcImageS, {bands: ['pc1', 'pc2', 'pc6'], min: -2, max: 2}, 'Sentinel 2 - PCA  126');
*/

/*var g =ee.Image.cat([pcImageL.select('pc2').multiply(1.50),pcImageS.select('pc2').multiply(0.375),pcImageS.select('pc5').multiply(1.125)])
Map.addLayer(g,{min: -3, max: 3},'FCCs combination 3')
*/
var h =ee.Image.cat([pcImageLDS.select('pc2').multiply(1),pcImageADS.select('pc3').multiply(1),pcImageS.select('pc5').multiply(1)])
//Map.addLayer(h,{min: -3, max: 3},'FCCs combination 1')
var h =ee.Image.cat([pcImageLDS.select('pc2').multiply(1.0),pcImageS.select('pc2').multiply(1.0),pcImageA.select('pc2').multiply(1)])
Map.addLayer(h,{min: -3, max: 3},'FCCs combination 2')

//SDS-1,5, ADS-3, LDS-2,3,5, A-2, S-2,5, L-1,3
var h =ee.Image.cat([pcImageSDS.select('pc5').multiply(1),pcImageLDS.select('pc5').multiply(1.0),pcImageA.select('pc2').multiply(1)])
Map.addLayer(h,{min: -3, max: 3},'FCCs combination 3')

var i =ee.Image.cat([pcImageS.select('pc2').multiply(1),pcImageL.select('pc1').multiply(1.5),pcImageSDS.select('pc5').multiply(0.5)])
Map.addLayer(i,{min: -3, max: 3},'FCCs combination 4')

var i =ee.Image.cat([pcImageS.select('pc2').multiply(1),pcImageA.select('pc2').multiply(1),pcImageL.select('pc3').multiply(1)])
Map.addLayer(i,{min: -3, max: 3},'FCCs combination 6')
/*Export.image.toDrive({
  image: g,
  description: 'FinalMap',
  folder: "TestExports",
  scale: 15,
  region: region
});
*/


Map.setCenter(72.13059794327728,34.550499453749666, 12);



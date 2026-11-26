// St John / Mt Finlay --- independent GFD flood observation discovery and export.
// Run in https://code.earthengine.google.com/ . No Hexalink simulation.
// Based on existing St John screening domain, simplified 700 m only for selection/export.
// GFD 2000–2018 mapped EVENTS, MODIS-derived, high cloud uncertainty.
// Reference: Tellman et al., Nature 2021, doi:10.1038/s41586-021-03695-w
// Bands: flooded 0/1, clear_views, jrc_perm_water. Export: 0 dry / 1 flood / 255 invalid.
// If no qualifying events, it prints that fact and does NOT create an export.

var aoi = ee.Geometry.MultiPolygon([[[[-9.838143,6.393358],[-9.839048,6.394261],[-9.838144,6.394263],[-9.838143,6.393358]]],[[[-9.83543,6.393363],[-9.828213,6.403324],[-9.83185,6.415982],[-9.843609,6.418676],[-9.834577,6.425022],[-9.838209,6.434062],[-9.828279,6.444028],[-9.811982,6.43139],[-9.785783,6.446808],[-9.77225,6.467633],[-9.782204,6.472141],[-9.772284,6.490247],[-9.788571,6.495649],[-9.790403,6.510118],[-9.81573,6.512792],[-9.822088,6.529967],[-9.839285,6.53808],[-9.860076,6.531713],[-9.858305,6.553424],[-9.87369,6.558825],[-9.875525,6.574198],[-9.89,6.575982],[-9.906366,6.622986],[-9.925391,6.638328],[-9.904636,6.666406],[-9.90287,6.689926],[-9.886607,6.701714],[-9.89115,6.71256],[-9.900194,6.710734],[-9.891183,6.73065],[-9.903847,6.728817],[-9.908401,6.74509],[-9.949169,6.772147],[-9.951913,6.787518],[-9.979066,6.791081],[-9.976386,6.808272],[-9.968244,6.809193],[-9.97643,6.829979],[-9.924913,6.863547],[-9.932185,6.879813],[-9.9159,6.882558],[-9.901458,6.903389],[-9.88786,6.891656],[-9.868874,6.902545],[-9.852558,6.889007],[-9.84533,6.896256],[-9.83265,6.890851],[-9.832683,6.909846],[-9.819137,6.927959],[-9.809173,6.923454],[-9.774853,6.969641],[-9.726904,6.988713],[-9.722399,7.002287],[-9.709711,6.993261],[-9.694326,6.996902],[-9.699794,7.021317],[-9.67446,7.032208],[-9.647292,7.027724],[-9.65186,7.056663],[-9.626529,7.072979],[-9.637415,7.087437],[-9.620218,7.092888],[-9.622043,7.10374],[-9.598532,7.129098],[-9.526104,7.145468],[-9.505328,7.193433],[-9.49264,7.18621],[-9.460059,7.210667],[-9.448278,7.204347],[-9.431984,7.214312],[-9.390306,7.200781],[-9.380354,7.213454],[-9.358614,7.212567],[-9.359532,7.227944],[-9.34324,7.244239],[-9.311535,7.243357],[-9.313368,7.275016],[-9.291625,7.273221],[-9.292539,7.28498],[-9.269897,7.294944],[-9.289844,7.321165],[-9.273552,7.347407],[-9.283516,7.342878],[-9.286241,7.353732],[-9.265407,7.361885],[-9.259072,7.375458],[-9.220113,7.38362],[-9.209249,7.403526],[-9.190219,7.404439],[-9.187509,7.422531],[-9.174823,7.42706],[-9.177547,7.440627],[-9.15943,7.462344],[-9.152179,7.459633],[-9.142217,7.48225],[-9.129529,7.484063],[-9.125001,7.498538],[-9.13316,7.503058],[-9.097816,7.525682],[-9.090571,7.555534],[-9.065191,7.570917],[-9.048877,7.60077],[-9.067014,7.633332],[-9.052514,7.682181],[-9.042539,7.676754],[-9.036192,7.686705],[-9.020777,7.685802],[-9.015337,7.700275],[-9.006269,7.69123],[-8.964557,7.69756],[-8.94189,7.678562],[-8.929192,7.694842],[-8.890203,7.686693],[-8.869342,7.703874],[-8.835791,7.702053],[-8.82854,7.693004],[-8.772297,7.735493],[-8.717875,7.751742],[-8.707872,7.790631],[-8.669756,7.819549],[-8.671561,7.831309],[-8.644335,7.847569],[-8.613453,7.891864],[-8.565386,7.878249],[-8.537261,7.880933],[-8.530941,7.854694],[-8.516438,7.845632],[-8.500096,7.857372],[-8.465625,7.858234],[-8.451144,7.833794],[-8.452084,7.808468],[-8.433973,7.786736],[-8.464822,7.776826],[-8.470277,7.765979],[-8.497487,7.763298],[-8.51021,7.7407],[-8.502968,7.729837],[-8.517499,7.710858],[-8.510254,7.702709],[-8.539333,7.644851],[-8.49588,7.587816],[-8.430551,7.632962],[-8.40881,7.62027],[-8.517712,7.52],[-8.530433,7.491069],[-8.503298,7.442193],[-8.508751,7.428631],[-8.484297,7.416844],[-8.492482,7.390622],[-8.506991,7.380688],[-8.503396,7.354452],[-8.523322,7.360806],[-8.537841,7.340016],[-8.520641,7.32643],[-8.507087,7.293852],[-8.515254,7.281197],[-8.527928,7.289351],[-8.543344,7.273989],[-8.549676,7.283041],[-8.564177,7.276724],[-8.566919,7.251399],[-8.583237,7.237845],[-8.575096,7.225174],[-8.62585,7.18994],[-8.616803,7.177269],[-8.624981,7.144711],[-8.649444,7.131161],[-8.644938,7.102212],[-8.613288,7.048817],[-8.619646,7.023494],[-8.599731,7.019859],[-8.575305,6.998128],[-8.607033,6.948403],[-8.63239,6.936664],[-8.643277,6.905012],[-8.653235,6.902306],[-8.647825,6.87426],[-8.684944,6.860717],[-8.699413,6.881532],[-8.708474,6.86435],[-8.727473,6.879739],[-8.755534,6.876136],[-8.750117,6.8499],[-8.758269,6.83724],[-8.748327,6.811001],[-8.761026,6.753114],[-8.750171,6.745871],[-8.760131,6.732307],[-8.732095,6.70606],[-8.741152,6.688877],[-8.753818,6.691597],[-8.759259,6.665367],[-8.74118,6.635506],[-8.745712,6.619225],[-8.706811,6.621013],[-8.718581,6.603832],[-8.70864,6.586639],[-8.722219,6.568555],[-8.740318,6.554091],[-8.740309,6.571278],[-8.763827,6.572195],[-8.805455,6.527887],[-8.83439,6.555035],[-8.851575,6.55685],[-8.856102,6.542378],[-8.881427,6.55143],[-8.893189,6.538768],[-8.924846,6.546011],[-8.921229,6.537869],[-8.961029,6.524304],[-8.964648,6.506212],[-8.991782,6.517974],[-9.028863,6.481788],[-9.091263,6.454642],[-9.085824,6.375942],[-9.131042,6.403071],[-9.140987,6.395831],[-9.155458,6.404873],[-9.166307,6.395824],[-9.181681,6.400341],[-9.191621,6.38315],[-9.186191,6.371392],[-9.209696,6.359623],[-9.235016,6.36504],[-9.250397,6.38403],[-9.258533,6.379502],[-9.250382,6.353273],[-9.290168,6.355966],[-9.288354,6.346921],[-9.297395,6.346011],[-9.32273,6.374943],[-9.334481,6.369508],[-9.343531,6.381262],[-9.363418,6.373107],[-9.367949,6.386673],[-9.386035,6.388469],[-9.400514,6.402931],[-9.403201,6.370364],[-9.421279,6.363112],[-9.425788,6.348635],[-9.443882,6.361284],[-9.466453,6.324176],[-9.488182,6.356721],[-9.489979,6.344055],[-9.508963,6.341323],[-9.52253,6.345833],[-9.52978,6.363013],[-9.548767,6.362088],[-9.556881,6.34037],[-9.574961,6.337636],[-9.585823,6.349383],[-9.596658,6.336707],[-9.625588,6.334864],[-9.624635,6.294159],[-9.64814,6.294131],[-9.645406,6.276947],[-9.658952,6.265171],[-9.693312,6.272363],[-9.694193,6.255176],[-9.705951,6.259683],[-9.711336,6.23073],[-9.726698,6.227995],[-9.72668,6.214427],[-9.737531,6.21803],[-9.748337,6.189069],[-9.780846,6.169122],[-9.804294,6.13562],[-9.832291,6.123818],[-9.840385,6.098478],[-9.869276,6.082151],[-9.898203,6.08934],[-9.923559,6.120955],[-9.899183,6.134564],[-9.902817,6.146317],[-9.879356,6.168064],[-9.893857,6.192462],[-9.875791,6.198824],[-9.892104,6.225028],[-9.886699,6.235891],[-9.906613,6.252139],[-9.896718,6.280195],[-9.910326,6.308212],[-9.884147,6.329965],[-9.895915,6.33899],[-9.89233,6.357087],[-9.874245,6.355308],[-9.857103,6.377046],[-9.859832,6.386991],[-9.83543,6.393363]],[[-8.92575,6.545107],[-8.926655,6.546011],[-8.926655,6.545107],[-8.92575,6.545107]],[[-8.738428,6.707872],[-8.737523,6.706967],[-8.737523,6.707872],[-8.738428,6.707872]]]], null, false);
Map.centerObject(aoi, 8);
Map.addLayer(aoi, {color:'orange'}, 'St John modeled-area approximation');

var gfd = ee.ImageCollection('GLOBAL_FLOOD_DB/MODIS_EVENTS/V1')
  .filterBounds(aoi);
print('GFD events intersecting the model AOI footprint:', gfd.size());

var stage = gfd.map(function(im) {
  var flooded = im.select('flooded');
  var hasClearObservation = im.select('clear_views').gt(0);
  var notPermanent = im.select('jrc_perm_water').eq(0);
  var valid = flooded.mask().and(hasClearObservation).and(notPermanent);
  var floodedValid = flooded.updateMask(valid);
  // Count at 250 m, the effective MODIS observation scale.
  var stats = floodedValid.reduceRegion({
    reducer: ee.Reducer.sum().combine({reducer2:ee.Reducer.count(),sharedInputs:true}),
    geometry:aoi,scale:250,crs:'EPSG:32629',maxPixels:10000000,tileScale:4
  });
  return im.set({
    'roi_flood_pixels': ee.Number(ee.Algorithms.If(stats.get('flooded_sum'), stats.get('flooded_sum'), 0)),
    'roi_valid_pixels': ee.Number(ee.Algorithms.If(stats.get('flooded_count'), stats.get('flooded_count'), 0))
  });
});
// A minimum of 10 observed flooded pixels *in this AOI* avoids export of empty events.
var eligible = stage.filter(ee.Filter.gte('roi_flood_pixels',10))
  .filter(ee.Filter.gte('roi_valid_pixels',100));
var listings = ee.FeatureCollection(eligible.toList(100).map(function(item) {
  var im=ee.Image(item);
  return ee.Feature(null, {
    'id':im.get('id'),
    'start_date':ee.Date(im.get('system:time_start')).format('YYYY-MM-dd'),
    'roi_flood_pixels':im.get('roi_flood_pixels'),
    'roi_valid_pixels':im.get('roi_valid_pixels'),
    'countries':im.get('countries')
  });
}));
print('Eligible event inventory (review before using as validation):',listings);

eligible.size().evaluate(function(n) {
  print('QUALIFYING EVENTS IN ST JOHN AOI:', n);
  if (!n) {
    print('No GFD event with qualifying observed flood in AOI. NO EXPORT made. Use FLOD-Africa/Copernicus archive instead.');
    return;
  }
  // Predefined selection: event with most observed flood pixels inside model AOI.
  // This selects observational support, NOT a score-optimizing Hexalink configuration.
  var im = ee.Image(eligible.sort('roi_flood_pixels', false).first());
  print('Selected GFD source event metadata:', im.toDictionary([
    'id','system:time_start','roi_flood_pixels','roi_valid_pixels','countries'
  ]));
  // Event metadata exported alongside the corresponding 250-m flood raster.
  // Both tasks use the same selected ee.Image, preserving provenance.
  var eventStart = ee.Date(im.get('system:time_start')).format('YYYY-MM-dd');
  var metadata = ee.FeatureCollection([ee.Feature(null, {
    'collection': 'GLOBAL_FLOOD_DB/MODIS_EVENTS/V1',
    'event_id': im.get('id'),
    'system_index': im.get('system:index'),
    'event_start_date_UTC': eventStart,
    'event_start_time_ms': im.get('system:time_start'),
    'roi_flood_pixels': im.get('roi_flood_pixels'),
    'roi_valid_pixels': im.get('roi_valid_pixels'),
    'selection': 'most independently observed flooded pixels in St John AOI',
    'corresponding_raster': 'StJohn_GFD_Observed_Event_250m.tif'
  })]);
  Export.table.toDrive({
    collection: metadata,
    description: 'StJohn_GFD_Event_Metadata',
    folder: 'Hexalink_StJohn_Flood_Validation',
    fileNamePrefix: 'StJohn_GFD_Event_Metadata',
    fileFormat: 'CSV',
    selectors: ['collection','event_id','system_index','event_start_date_UTC',
      'event_start_time_ms','roi_flood_pixels','roi_valid_pixels','selection',
      'corresponding_raster']
  });
  print('METADATA CSV TASK CREATED: start both CSV and GeoTIFF tasks in Tasks.');
  var flooded=im.select('flooded');
  var valid=flooded.mask()
    .and(im.select('clear_views').gt(0))
    .and(im.select('jrc_perm_water').eq(0));
  var observed=flooded.updateMask(valid).clip(aoi).unmask(255).toUint8().rename('observed_flood');
  Map.addLayer(observed.updateMask(observed.eq(1)),{palette:['0088ff']},'Satellite-observed event flood');
  Export.image.toDrive({
    image:observed,
    description:'StJohn_GFD_Observed_Event_250m',
    folder:'Hexalink_StJohn_Flood_Validation',
    fileNamePrefix:'StJohn_GFD_Observed_Event_250m',
    region:aoi,
    crs:'EPSG:32629',
    scale:250,
    maxPixels:10000000,
    fileFormat:'GeoTIFF',
    formatOptions:{cloudOptimized:true}
  });
  print('ONE GeoTIFF export task created. Start it in Tasks; also record printed event ID and start date.');
});

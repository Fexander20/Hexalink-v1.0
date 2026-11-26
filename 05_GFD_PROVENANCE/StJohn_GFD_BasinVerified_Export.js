// HEXALINK — exact St. John GFD event provenance
// This reproduces the ORIGINAL basin-specific AOI and event-selection rules.
// The AOI vertices are losslessly encoded to six decimal places for easy copying.
// Creates a matched CSV + GeoTIFF. No model results are changed.

var basinRings = [["|bnwQ{_feKpw@mw@ow@CApw@"], ["jyhwQe`feKabMqmRhbFcvW|}UkgDosPskK~aF_tPskR{mRqy^ztWmdr@sb]ylYqtg@bmRwxG_kRsjb@|x^spInqBig[|mp@cfDjlKmp`@xq`@azNlrg@|lKumB}ki@p`]qpItqBy_]tg[onBz}^wxzA`dd@{}\\epg@{yu@kmB_}l@mw^w_V|zG{dTftPbqBerPw{e@nvWpqBr{Gax^~rnAazs@njDu_]``t@u}EofDmq`@{{Nqx@r~Ncrg@yrcB_q`AneMsw^yx^qjDse[}tg@{pYh|Usad@qgTwz^bmYwbMadMowWxpI`Aebd@smYakb@wmRpxG_`bAueyAys|A_gd@qxGkoY_xWbsPq`]qbFvtI}tn@knp@ugT_at@fwGn|Guow@enp@wz^jgTsf[yq`@usI`qBgeTm|l@{op@wmlCc~^oqg@yt|A_xWlbMis~@qwn@i_V~iKky^ymR{kpAtlY_mRawWwmi@lv@jx@a`]gy^my^q||@bv@pqBuy|@}mi@doBbx@}}Ucfk@wmRt}e@yer@gy^cgr@vmR`zGhiDkeTcug@q|N}jKioY}akAc}N_fTc{e@kdd@ax@khDwib@{wWazGfiD}nYikb@ili@edMlhDsmRqdk@_xWipB_zGug[||NoyG_`dA_ek@ycMwhy@gqp@m`]sz^yhy@plb@cr~@gi[al~AmnRdrIukK}lRmb]lw@_sIqg[wuPhtP_npAsjKugk@jbd@sxWox^yckAh|Nyvg@yp`@}o`AxpBedMptPezlBq~qAkhiBqv^epRq}jAgmiAknw@xoB_~Usdt@gw^ci{@mouAe{|A|qYy|u@wfD_jK|fr@mi[juPk|^w|UmibA{t@ah[nvn@vy@zmp@}jb@fmi@`g{@jjR|sI|dTrct@pfDdzWjck@scM|eTdk[dad@ycMh|Nlxw@b_pByzsAtknBar~BsdwAymi@fxWjusEzybE`zWdow@}~s@vm~AxsIrnYkwn@t_Vp~Nzer@xi[zkRu_Fvfr@j|e@clKlj[jrg@_r`@bpYcnYbs~@l}N|uWbwWs|Nnb]b_]vjKwtPhi[xiKjjDxmp@zz^bnYy{N|vWbcbBbycAmtP|vWb~Nzq~@|wn@zmYsxGhpw@cy|@dhgBjlKtmp@u{e@dbFsun@dmi@~}|@xb`Bxop@t|UlgTfy|@jmRbhDcqIzwu@|ngAlmYhg[}sg@huPzp`@lbd@y`]xxu@d`FqqIvfr@n|NfvWklR|fr@txW|`pBmeTtcMnmRvnYgwu@lgr@`uP|p`@rvW_iD`sIjfr@}hb@hiy@fzGpx^i~jAwnBr~Uxp`@ilRpq`@toYfib@djb@~f[Qeq`@z|l@ix@vhpAfpuAlow@w_t@`q`@mpB|yGng[xmp@wtPb~UjvWpy|@ucMaaFz{NnvlAxnYdaFvib@z~s@c~UplgArteA~zxBr_t@}rIvuxCbiwAq~s@plRncMlg[ctP`eTptPz_]iyGflRlq`@krIz}U`|l@p~Unmp@qqIh`]{ad@n{N~yGm|Nhfr@rulAigDkpBhtP`tPzw@lnp@gow@l}UtrIrtPs}U|ye@t|NdzG{nYjib@goB|g[{f[|fDlr~@zhb@fdMxxGxg[zib@quWtak@fngA`mi@aq~@hoBrvWnad@viD|nY{xGbdMwp`@tad@xx@bzNjli@~hb@ziDzeTe}UddTfwWbow@drBqz@`onA`|l@v@{iD~p`@rmY~~UnbbAo`M`v@dq`@z}UuxGpoIppw@b_]|iDc@~nYdeTe`FjbT`qw@xn~@t}e@nxl@zl`Axtu@r`VzxNvnp@tlw@l{^|nw@i`Mvop@}v|@orn@qqYbbFq}Uiyl@eni@hi[{sn@chb@slKpz^wdr@ypI}eTr{e@ov^miRoxu@nqYavu@ecr@qni@n~UasPa_Fajb@iib@dnBkn`@smi@piDqlRctn@gmK", "jbx_PelnnKpw@ow@?nw@qw@?", "vfjtO_ilxKqw@pw@?qw@pw@?"]];

function decodeRing(s) {
  var i = 0, x = 0, y = 0, points = [];
  function nextDelta() {
    var shift = 0, value = 0, b;
    do {
      b = s.charCodeAt(i++) - 63;
      value |= (b & 31) << shift;
      shift += 5;
    } while (b >= 32);
    return (value & 1) ? -((value >> 1) + 1) : (value >> 1);
  }
  while (i < s.length) {
    x += nextDelta();
    y += nextDelta();
    points.push([x / 1000000, y / 1000000]);
  }
  return points;
}

var coords = basinRings.map(function(poly) {
  return poly.map(function(ring) { return decodeRing(ring); });
});
var aoi = ee.Geometry.MultiPolygon(coords, null, false);
Map.centerObject(aoi, 8);
Map.addLayer(aoi, {color: 'orange'}, 'Original St. John AOI');

var gfd = ee.ImageCollection('GLOBAL_FLOOD_DB/MODIS_EVENTS/V1')
  .filterBounds(aoi);

var assessed = gfd.map(function(im) {
  var flood = im.select('flooded');
  var valid = flood.mask()
    .and(im.select('clear_views').gt(0))
    .and(im.select('jrc_perm_water').eq(0));
  var stats = flood.updateMask(valid).reduceRegion({
    reducer: ee.Reducer.sum().combine({
      reducer2: ee.Reducer.count(), sharedInputs: true
    }),
    geometry: aoi,
    scale: 250,
    crs: 'EPSG:32629',
    maxPixels: 10000000,
    tileScale: 4
  });
  return im.set({
    roi_flood_pixels: ee.Number(ee.Algorithms.If(
      stats.get('flooded_sum'), stats.get('flooded_sum'), 0)),
    roi_valid_pixels: ee.Number(ee.Algorithms.If(
      stats.get('flooded_count'), stats.get('flooded_count'), 0))
  });
});

var eligible = assessed
  .filter(ee.Filter.gte('roi_flood_pixels', 10))
  .filter(ee.Filter.gte('roi_valid_pixels', 100));
print('Eligible events in St. John AOI:', eligible.size());

eligible.size().evaluate(function(n) {
  if (!n) {
    print('No qualifying GFD events. Nothing exported.');
    return;
  }
  // ORIGINAL criterion: maximum observed flood within the basin AOI.
  var selected = ee.Image(eligible.sort('roi_flood_pixels', false).first());
  print('SELECTED EVENT — PLEASE COPY:', selected.toDictionary([
    'id', 'system:index', 'system:time_start',
    'roi_flood_pixels', 'roi_valid_pixels', 'countries'
  ]));
  print('Does it match DFO_3982 from your earlier CSV?',
    ee.Number(selected.get('id')).eq(3982));

  var metadata = ee.FeatureCollection([ee.Feature(null, {
    collection: 'GLOBAL_FLOOD_DB/MODIS_EVENTS/V1',
    event_id: selected.get('id'),
    event_index: selected.get('system:index'),
    event_start_date: ee.Date(selected.get('system:time_start')).format('YYYY-MM-dd'),
    roi_flood_pixels: selected.get('roi_flood_pixels'),
    roi_valid_pixels: selected.get('roi_valid_pixels'),
    event_id_is_DFO_3982: ee.Number(selected.get('id')).eq(3982),
    region: 'Original St John basin-specific AOI',
    paired_raster: 'StJohn_GFD_BasinVerified_250m.tif'
  })]);

  var flood = selected.select('flooded');
  var valid = flood.mask()
    .and(selected.select('clear_views').gt(0))
    .and(selected.select('jrc_perm_water').eq(0));
  var observed = flood.updateMask(valid)
    .clip(aoi).unmask(255).toUint8().rename('observed_flood');
  Map.addLayer(observed.updateMask(observed.eq(1)),
    {min: 0, max: 1, palette: ['0088ff']}, 'Observed flood');

  Export.table.toDrive({
    collection: metadata,
    description: 'StJohn_GFD_BasinVerified_Metadata',
    folder: 'Hexalink_StJohn_Flood_Validation',
    fileNamePrefix: 'StJohn_GFD_BasinVerified_Metadata',
    fileFormat: 'CSV'
  });

  Export.image.toDrive({
    image: observed,
    description: 'StJohn_GFD_BasinVerified_250m',
    folder: 'Hexalink_StJohn_Flood_Validation',
    fileNamePrefix: 'StJohn_GFD_BasinVerified_250m',
    region: aoi,
    crs: 'EPSG:32629',
    scale: 250,
    maxPixels: 10000000,
    fileFormat: 'GeoTIFF',
    formatOptions: {cloudOptimized: true}
  });
  print('Two tasks created: metadata CSV and matching flood GeoTIFF.');
});

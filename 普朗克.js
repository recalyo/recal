var geometry = ee.Geometry.MultiPoint();
var roi = ee.FeatureCollection("projects/downloadgee-454605/assets/CD_wuchengqu_dissolve");

// 1. 云掩膜函数
// 利用Landsat 8获取地表温度 // 定义一个函数用于云和阴影的掩膜处理
function maskL8sr(image) {
   // 创建QA_PIXEL质量掩膜
   var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0); 
   // 饱和度掩膜，去除云
   var saturationMask = image.select('QA_RADSAT').eq(0);
   // 对光学波段进行归一化，计算正常值
   var opticalBands = image.select('SR_B.*').multiply(0.0000275).add(-0.2);
   // 对热红外波段进行归一化，获取地表温度映射
   var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
   // 将处理后的波段添加到影像中，并更新掩膜
   return image.addBands(opticalBands, null, true)
               .addBands(thermalBands, null, true)
               .updateMask(qaMask)
               .updateMask(saturationMask);}
               
// 加载Landsat 8 Collection 2 Level-2 Surface Reflectance数据
var landsat8_sr = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                     .filterBounds(roi)  // 选择感兴趣区域
                     .map(maskL8sr) // 处理影像以去除云和阴影 
                     .filterDate('2023-07-01', '2023-08-31') // 选择时间范围
                     .median()  // 取该区域该时间范围内的中值影像
                     .clip(roi);
 
// 选择影像
var image2 = landsat8_sr;  
 
// 获取红色波段和近红外波段
var red = image2.select('SR_B4');  // 红色波段 (Band 4)
var nir = image2.select('SR_B5');  // 近红外波段 (Band 5)
var trad = image2.select('ST_TRAD').multiply(0.001);  // 转换为辐射温度
var urad = image2.select('ST_URAD').multiply(0.001);  // 转换为大气辐射
var atran = image2.select('ST_ATRAN').multiply(0.0001);  // 大气透过率
var drad = image2.select('ST_DRAD').multiply(0.001);  // 辐射强度
 
// 计算NDVI
var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
 
// 计算地表比辐射率（FVC），基于NDVI进行分类计算
/*
var fvc = ndvi.lte(0).multiply(0.995)
    .add(ndvi.gt(0).and(ndvi.lt(0.7))
        .multiply(0.9589 + 0.086 * nir.subtract(red)
        .subtract(0.0671 * nir.subtract(red).multiply(nir.subtract(red))))
    )
    .add(ndvi.gte(0.7)
        .multiply(0.9625 + 0.0614 * nir.subtract(red)
        .subtract(0.0461 * nir.subtract(red).multiply(nir.subtract(red))))
    );
*/
var fvc = ee.Image().expression(
'(b1 > 0.7) * 1 + (b1 < 0.05) * 0 + (b1 >= 0.05 && b1 <= 0.7) * ((b1 - 0.05) / (0.7 - 0.05))', 
  {
    b1: ndvi 
  });
 
// 计算黑体辐射（BHR）
var bhr = trad.subtract(urad).subtract(atran.multiply(trad.subtract(urad).multiply(1.0 - fvc).multiply(drad)))
    .divide(atran.multiply(fvc));  // 计算黑体辐射
 /*   
var bhr = ee.Image().expression(
  '(b1 <= 0) * 0.995 + ' +
  '(b1 > 0 && b1 < 0.7) * (0.9589 + 0.086 * b2 - 0.0671 * b2 * b2) + ' +
  '(b1 >= 0.7) * (0.9625 + 0.0614 * b2 - 0.0461 * b2 * b2)', 
  {
    b1: ndvi,  
    b2: fvc  
  });*/
 
// 计算地表温度（LST）
var lst_Lc = ee.Image(1321.08)
    .divide((ee.Image(774.89).divide(trad)).add(1).log())
    .subtract(273.15);  // 计算LST，单位为摄氏度
 
 
 
 
// 将处理后的LST图层添加到地图上
Map.addLayer(lst_Lc, { 
    'min': 2,  
    'max': 49,  
     'palette': [ 
    "eff3ff", "c6dbef", "9ecae1", "6baed6",
    "4292c6", "2171b5", "084594", "fff5f0",
    "fee0d2", "fcbba1", "fc9272", "fb6a4a", 
    "ef3b2c", "cb181d", "99000d" 
    ]}, 'lst_Lc');
    
    
print("产品处理后直方图", ui.Chart.image.histogram(lst_Lc, roi, 100, 258));
print(lst_Lc);
 
function exportImage(image, roi, fileName) {
  Export.image.toDrive({
    image: image,
    description: "Landsat8_svm" + fileName,   // 导出描述
    fileNamePrefix: fileName,  // 文件前缀
    folder: "Landsat8",  // 保存的文件夹 
    scale: 30,  // 分辨率
    region: roi,  // 导出区域
    maxPixels: 1e13,  // 最大像素数量 
    crs: "EPSG:4326"  // 设置投影
    });  }// 调用导出函数，将Landsat 8的地表温度影像导出到Google Drive
exportImage(lst_Lc, roi, "lst_Lc");
// Remove this if you don't need to generate new cell intensity measurements (it may be quite slow)
// 0.2.0
//Original source https://forum.image.sc/t/stardist-qupath-cell-segmentation-transfer-on-aligned-image-set/45362/5
//With substitutions from @petebankhead and @mesencephalon
// edited by Yulan weng

import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerMetadata
import qupath.lib.images.servers.TransformedServerBuilder

def imageData = getCurrentImageData()
def server = new TransformedServerBuilder(imageData.getServer())
    .deconvolveStains(imageData.getColorDeconvolutionStains())
    .build()
 // or server = getCurrentServer() //for IF   
def measurements = ObjectMeasurements.Measurements.values() as List
def compartments = ObjectMeasurements.Compartments.values() as List // Won't mean much if they aren't cells...
def downsample = 1.0
for (detection in getDetectionObjects()) {
  ObjectMeasurements.addIntensityMeasurements(
      server, detection, downsample, measurements, compartments
      )
}

// add cell label
i = 1
getCellObjects().each{
it.getMeasurementList().putMeasurement("Label", i)
i++
}
println("done")

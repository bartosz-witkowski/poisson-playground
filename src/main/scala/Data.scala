package poisson

import javax.imageio.ImageIO

import java.awt.image.{DataBufferByte, DataBufferInt, BufferedImage}
import scala.concurrent.{Future, Await}
import scala.concurrent.duration.Duration
import scala.concurrent.ExecutionContext.Implicits.global

import java.io.File

class RawImageData(
  val width: Int,
  val height: Int,
  val r: Array[Double],
  val g: Array[Double],
  val b: Array[Double]) {

  def transformWH(f: (Int, Int, Array[Double]) => Array[Double]): RawImageData = {
    val rFut = Future(f(width, height, r))
    val gFut = Future(f(width, height, g))
    val bFut = Future(f(width, height, b))

    val rr = Await.result(rFut, Duration.Inf)
    val gg = Await.result(gFut, Duration.Inf)
    val bb = Await.result(bFut, Duration.Inf)

    RawImageData(width, height, rr, gg, bb)
  }

  def transform(f: Array[Double] => Array[Double]): RawImageData = {
    transformWH { (_, _, c) =>
      f(c)
    }
  }
}

object RawImageData {
  /*
   * original:  padded:
   * ab         0000
   * cd         0ab0
   *            0cd0
   *            0000
   *
   */
  def fromBufferedImagePadded(img: BufferedImage): RawImageData = {
    from(img, 1)
  }

  def fromBufferedImage(img: BufferedImage): RawImageData = {
    from(img, 0)
  }

  private def from(
      img: BufferedImage,
      padding: Int): RawImageData = {
    val pixels = (img.getRaster().getDataBuffer().asInstanceOf[DataBufferByte]).getData()
    
    val height = img.getHeight + (padding * 2)
    val width = img.getWidth + (padding * 2)

    val size = width * height
    val b = new Array[Double](size)
    val g = new Array[Double](size)
    val r = new Array[Double](size)
    val (stride, rOff, gOff, bOff) = if (img.getAlphaRaster() == null) {
      (3, 2, 1, 0) 
    } else {
      (4, 3, 2, 1)
    }

    var src = 0
   
    var y = padding
    while (y < height - padding) {
      
      var x = padding
      while (x < width - padding) {
        val dest = x + (y * width)
        
        b(dest) = java.lang.Byte.toUnsignedInt(pixels(src + bOff)) / 255.0
        g(dest) = java.lang.Byte.toUnsignedInt(pixels(src + gOff)) / 255.0
        r(dest) = java.lang.Byte.toUnsignedInt(pixels(src + rOff)) / 255.0
        
        src += stride

        x += 1
      }

      y += 1
    }

    RawImageData(width, height, r, g, b)
  }
}

def writeImageToPngFile(img: BufferedImage, file: File): Unit = {
  ImageIO.write(img, "png", file)
}

def writeImageToPngFile(img: RawImageData, file: File): Unit = {
  writeImageToPngFile(toBufferedImage(img), file)
}

def toBufferedImage(img: RawImageData): BufferedImage = {
  val channs = Array(img.b, img.g, img.r)

  inline def get(offset: Int, chann: Int) = {
    val d: Double = channs(chann)(offset)

    val i = (Math.round(255.0 * d)).toInt
    if (i < 0) 0 else 
    if (i > 255) 255 else
    i
  }

  var offset = (img.width * img.height)
  val source = new Array[Int](offset)
  offset -= 1

  while (offset >= 0) {
    source(offset) = (
       get(offset, 0) |
      (get(offset, 1) << 8) | // r
      (get(offset, 2) << 16)
    )
    offset -= 1
  }

  val bufferedImage = new BufferedImage(img.width, img.height, BufferedImage.TYPE_INT_RGB)

  val dest = bufferedImage.getRaster().getDataBuffer.asInstanceOf[DataBufferInt].getData()
  System.arraycopy(source, 0, dest, 0, source.length)

  bufferedImage
}


case class ImageLaplacian(img: RawImageData)

/*
 * Single channel data to grayscale image with normalization
 */
def dataToImage(
  b: Array[Double],
  width: Int,
  height: Int): BufferedImage = {

  val min = b.min
  val max = b.max
  val span = max - min

  inline def toInt(unnormalized: Double) = {
    val d = (unnormalized - min) / span

    val i = (Math.round(255.0 * d)).toInt
    if (i < 0) 0 else 
    if (i > 255) 255 else
    i
  }

  var offset = (width * height)
  val source = new Array[Int](offset)
  offset -= 1

  while (offset >= 0) {
    source(offset) = (
       toInt(b(offset)) |
      (toInt(b(offset)) << 8) | 
      (toInt(b(offset)) << 16)
    )

    offset -= 1
  }

  val img = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)

  val dest = img.getRaster().getDataBuffer.asInstanceOf[DataBufferInt].getData()
  System.arraycopy(source, 0, dest, 0, source.length)

  img
}

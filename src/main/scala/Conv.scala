package poisson

case class Kernel(size: Int, data: Array[Double]) {
  val offset = ((size - 1) / 2)

  def :*(const: Double): Kernel = {
    Kernel(size, data.map(_ * const))
  }

  // dx goes from   ... -1, 0, 1 ...
  inline def fromDeltas(dx: Int, dy: Int): Double = {
    data((dx + offset) + ((dy + offset) * size))
  }
}
object Kernel {
  def apply(xs: Double*): Kernel = {
    val size  = math.sqrt(xs.size).toInt

    if ((size * size) != xs.size || size < 2) {
      throw new IllegalArgumentException("Kernel must have size 3, 5, 7, ... ")
    }

    Kernel(size, xs.toArray)
  }
}

/*
 * Performs one channel convolution on dataIn storing results in dataOut. 
 *
 * The behavior near the border treats all border-pixels as if they had infinite width. 
 */
def convExtended_!(
    width: Int, 
    height: Int, 
    dataIn: Array[Double], 
    dataOut: Array[Double], 
    kernel: Kernel): Unit = {

  val lim = (kernel.size - 1) / 2

  var y = 0
  while (y < height - 1) {
    var x = 0
    while (x < width - 1) {
      var sum = 0.0
      var ky = -lim

      while (ky <= lim) {
        var kx = -lim
        while (kx <= lim) {
          val sx = math.min(width - 1,  math.max(0, x + kx))
          val sy = math.min(height - 1, math.max(0, y + ky))

          sum += dataIn(sx + (sy * width)) * kernel.fromDeltas(kx, ky)
          kx += 1
        }

        ky += 1
      }

      dataOut(x + (y * width)) = sum

      x += 1
    }

    y += 1
  }
}

def convExtended(img: RawImageData, kernel: Kernel): RawImageData = {
  img.transform { chann => 
    convExtended(img.width, img.height, chann, kernel)
  }
}

def convExtended(
    width: Int, 
    height: Int, 
    dataIn: Array[Double], 
    kernel: Kernel): Array[Double] = {
  val output = new Array[Double](width * height)
  convExtended_!(width, height, dataIn, output, kernel)
  output
}


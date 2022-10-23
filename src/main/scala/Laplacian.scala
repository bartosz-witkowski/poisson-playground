package poisson

def fdx(width: Int, height: Int, buffer: Array[Double]): Array[Double] = {
  val result: Array[Double] = new Array[Double](width * height)
  
  var y = 0
  while (y < height) {
    var x = 1

    while (x < height - 1) {
      val offset = x + (y * width)
      result(offset) = buffer(offset) - buffer(offset - 1)

      x += 1
    }
    y += 1
  }

  result
}

def laplacianOf(width: Int, height: Int, buffer: Array[Double]): Array[Double] = {
  inline def b(x: Int, y: Int): Double = {
    buffer((y * width) + x)
  }

  val laplaceB: Array[Double] = new Array[Double](width * height)
  var dest = width + 1
  var y = 1

  while (y < height - 1) {
    var x = 1

    while (x < width - 1) {
      laplaceB(dest) = -b(x, y - 1) -b(x - 1, y) + (4*b(x, y)) - b(x + 1, y) -b(x, y + 1)
      dest += 1
      x += 1
    }
    dest += 2
    y += 1
  }

  laplaceB
}

def laplacianOf(rawImageData: RawImageData): ImageLaplacian = {
  import rawImageData.{width, height}

  ImageLaplacian(
    rawImageData.transform { c =>
      laplacianOf(rawImageData.width, rawImageData.height, c)
    }
  )
}

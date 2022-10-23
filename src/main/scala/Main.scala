package poisson

import javax.imageio.ImageIO
import java.io.File
import java.awt.image.{DataBufferByte, DataBufferInt, BufferedImage}

def getDiff(a: Array[Double], b: Array[Double]): Double = {
  var sum = 0.0
  a.indices.foreach { i =>
    sum = sum + math.abs(a(i) - b(i))
  }
  sum
}

def solve(lap: RawImageData, nIterations: Int = 50, nPre: Int = 20, nPost: Int = 60): RawImageData = {
  assert(nIterations / 2 > 0)
  import lap.{width, height}

  def solveChannel(chan: Array[Double]): Array[Double] = {
    val oldValues = new Array[Double](width * height)
    val newValues = new Array[Double](width * height)

    (0 until (nIterations / 2)).foreach { _ =>
      multigrid(width, height, oldValues, chan, newValues, nPre, nPost, 1.0)
      multigrid(width, height, newValues, chan, oldValues, nPre, nPost, 1.0)
    }

    oldValues
  }

  val t0 = System.currentTimeMillis
  val res = lap.transform(solveChannel)
  val t1 = System.currentTimeMillis
  println(s"Solved in ${t1 - t0}[ms]")

  res
}


def compare(): Unit =  {
  val img = ImageIO.read(new File("lenna.png"))

  val rawImageData = RawImageData.fromBufferedImagePadded(img)

  val laplacian: ImageLaplacian = laplacianOf(rawImageData)

  val oldB = rawImageData.b.clone()

  for {
    x <- 1 until rawImageData.width - 1
    y <- 1 until rawImageData.height - 1
  } oldB((y * rawImageData.width) + x) = 0.0

  val newB= oldB.clone()

  val t0 = System.currentTimeMillis
  (0 until 25).foreach { _ =>
    multigrid(rawImageData.width, rawImageData.height, oldB, laplacian.img.b, newB, 20, 60, 1.0)
    multigrid(rawImageData.width, rawImageData.height, newB, laplacian.img.b, oldB, 20, 60, 1.0)
  }
  val t1 = System.currentTimeMillis

  println(getDiff(newB, rawImageData.b) + " in " + (t1 - t0) + "[ms]")
}


def median(width: Int, height: Int, data: Array[Double]): Array[Double] = {
  inline def offset(x: Int, y: Int): Int = x + (y * width)
  val blured = new Array[Double](width * height)

  (1 until width - 1).foreach { x =>
    (1 until height -1).foreach { y =>
      val array = Array(
        data(offset(x - 1, y - 1)), data(offset(x, y - 1)), data(offset(x + 1, y - 1)),
        data(offset(x - 1, y)),     data(offset(x, y)),     + data(offset(x + 1, y)),
        data(offset(x - 1, y + 1)), data(offset(x, y + 1)), data(offset(x + 1, y + 1))
      ).sorted

      val median = array(array.size / 2)

      blured(offset(x, y)) = median
    }
  }

  blured
}

def normalize(width: Int, height: Int, data: Array[Double]): Array[Double] = {
  val min = data.min
  val max = data.max
  val span = max - min

  val out = new Array[Double](width * height)
  out.indices.foreach { i =>
    out(i) = (data(i) - min) / span
  }
  out
}

def nonlinearLap(l: ImageLaplacian)(f: (Int, Int, Array[Double]) => Array[Double]): RawImageData = {
  solve(l.img.transformWH(f)).transformWH(normalize)
}

def saveLaplacianOfFile(file: File): Unit = {
  val img = RawImageData.fromBufferedImagePadded(ImageIO.read(file))
  val laplacian: ImageLaplacian = laplacianOf(img)

  import img.{width, height}
  val imgLap = dataToImage(laplacian.img.b, width, height)
  val outName = file.getName.takeWhile(_ != '.') + "-laplacian.png"
  writeImageToPngFile(imgLap, new File(outName))
}

val BoxBlurKernel = Kernel(
  1, 1, 1,
  1, 1, 1,
  1, 1, 1
) :* (1/9.0)

val GaussianBlur = Kernel(
  1, 2, 1,
  2, 4, 2,
  1, 2, 1
) :* (1/16.0)

val Sharp = Kernel(
   0, -1,  0,
  -1,  5, -1,
   0, -1,  0
)

val Right = Kernel(
   0,  0,  0,
   0,  0,  1,
   0,  0,  0
)

val Emboss = Kernel(
  -2, -1,  0,
  -1,  1,  1,
   0,  1,  2)

def saveSymConvs(file: File): Unit = {
  val img = RawImageData.fromBufferedImagePadded(ImageIO.read(file))
  val laplacian: ImageLaplacian = laplacianOf(img)

  val name = file.getName.takeWhile(_ != '.')

  writeImageToPngFile(solve(convExtended(laplacian.img, BoxBlurKernel)), new File(s"${name}-lap-box-blur.png"))
  writeImageToPngFile(convExtended(img, BoxBlurKernel), new File(s"${name}-box-blur.png"))

  writeImageToPngFile(solve(convExtended(laplacian.img, GaussianBlur)), new File(s"${name}-lap-gauss-blur.png"))
  writeImageToPngFile(convExtended(img, GaussianBlur), new File(s"${name}-gauss-blur.png"))

  writeImageToPngFile(solve(convExtended(laplacian.img, Sharp)), new File(s"${name}-lap-sharp.png"))
  writeImageToPngFile(convExtended(img, Sharp), new File(s"${name}-sharp.png"))
}

def saveNonlinear(file: File): Unit = {
  val img = RawImageData.fromBufferedImagePadded(ImageIO.read(file))
  val laplacian: ImageLaplacian = laplacianOf(img)

  val name = file.getName.takeWhile(_ != '.')

  val rightSolved = solve(convExtended(laplacian.img, Right))
  writeImageToPngFile(rightSolved, new File(s"${name}-lap-right.png"))
  writeImageToPngFile(rightSolved.transformWH(normalize), new File(s"${name}-lap-right-norm.png"))
  writeImageToPngFile(convExtended(img, Right), new File(s"${name}-right.png"))

  val embossSolved = solve(convExtended(laplacian.img, Emboss))
  writeImageToPngFile(embossSolved, new File(s"${name}-lap-emboss.png"))
  writeImageToPngFile(embossSolved.transformWH(normalize), new File(s"${name}-lap-emboss-norm.png"))
  writeImageToPngFile(convExtended(img, Emboss), new File(s"${name}-emboss.png"))

  writeImageToPngFile(nonlinearLap(laplacian)(median), new File(s"${name}-lap-median.png"))
  writeImageToPngFile(img.transformWH(median), new File(s"${name}-median.png"))
}

def saveLowQuality(file: File): Unit = {
  val img = RawImageData.fromBufferedImagePadded(ImageIO.read(file))
  val laplacian: ImageLaplacian = laplacianOf(img)

  val name = file.getName.takeWhile(_ != '.')

  writeImageToPngFile(solve(laplacian.img, nIterations = 2, nPre = 2, nPost = 4).transformWH(normalize), new File(s"${name}-lq.png"))
}

def blendFaceOff(): Unit = {
  val target = RawImageData.fromBufferedImagePadded(ImageIO.read(Cage))
  val source = RawImageData.fromBufferedImage(ImageIO.read(TravoltaMouth))

  val blendX = 257
  val blendY = 239

  val sourceL: ImageLaplacian = laplacianOf(source)
  val targetL: ImageLaplacian = laplacianOf(target)

  def newChannel = new Array[Double](target.width * target.height)
  val r = newChannel
  val g = newChannel
  val b = newChannel
  val lr = newChannel
  val lg = newChannel
  val lb = newChannel
  
  val a = 0.5

  (0 until target.width - 1).foreach { x =>
    (0 until target.height - 1).foreach { y =>
      val offset = x + (y * target.width)

      @inline def inX = x >= (blendX + 1) && x < (blendX + source.width - 1)
      @inline def inY = y >= (blendY + 1) && y < (blendY + source.height - 1)
      @inline def inRect = inX && inY

      if (inRect) {
        val sourceOffset = (x - blendX) + ((y - blendY) * source.width)

        r(offset) = (a * source.r(sourceOffset))       + ((1 - a) * target.r(offset))
        g(offset) = (a * source.g(sourceOffset))       + ((1 - a) * target.g(offset))
        b(offset) = (a * source.b(sourceOffset))       + ((1 - a) * target.b(offset))

        lr(offset) = (a * sourceL.img.r(sourceOffset)) + ((1 - a) * targetL.img.r(offset))
        lg(offset) = (a * sourceL.img.g(sourceOffset)) + ((1 - a) * targetL.img.g(offset))
        lb(offset) = (a * sourceL.img.b(sourceOffset)) + ((1 - a) * targetL.img.b(offset))
      } else {
        r(offset) = target.r(offset)
        g(offset) = target.g(offset)
        b(offset) = target.b(offset)

        lr(offset) = targetL.img.r(offset)
        lg(offset) = targetL.img.g(offset)
        lb(offset) = targetL.img.b(offset)
      }
    }
  }

  val blended = RawImageData(target.width, target.height, r, g, b)
  val blendedL = RawImageData(target.width, target.height, lr, lg, lb)

  writeImageToPngFile(blended, new File(s"faceoff.png"))
  writeImageToPngFile(solve(blendedL), new File(s"faceoff-lap.png"))
}

def blendAll(f1: File, f2: File): Unit = {
  val alfa = RawImageData.fromBufferedImagePadded(ImageIO.read(f1))
  val beta = RawImageData.fromBufferedImagePadded(ImageIO.read(f2))

  if (alfa.width != beta.width) throw new RuntimeException(s"$f1.width must be $f2.width")
  if (alfa.height != beta.height) throw new RuntimeException(s"$f1.height must be $f2.height")
  
  import alfa.{width, height}

  val alfaL: ImageLaplacian = laplacianOf(alfa)
  val betaL: ImageLaplacian = laplacianOf(beta)

  def newChannel = new Array[Double](width * height)
  val r = newChannel
  val g = newChannel
  val b = newChannel

  val rr = newChannel
  val gg = newChannel
  val bb = newChannel
  
  val a = 0.3

  (0 until width - 1).foreach { x =>
    (0 until height - 1).foreach { y =>
      def maxabs(a: Double, b: Double): Double = if (a.abs > b.abs) { 
        a
      } else {
        b
      }

      val offset = x + (y * width)
      r(offset) = a * alfaL.img.r(offset) + (1 - a) * betaL.img.r(offset)
      g(offset) = a * alfaL.img.g(offset) + (1 - a) * betaL.img.g(offset)
      b(offset) = a * alfaL.img.b(offset) + (1 - a) * betaL.img.b(offset)

      rr(offset) = a * alfa.r(offset) + (1 - a) * beta.r(offset)
      gg(offset) = a * alfa.g(offset) + (1 - a) * beta.g(offset)
      bb(offset) = a * alfa.b(offset) + (1 - a) * beta.b(offset)
    }
  }

  val blendedL = RawImageData(width, height, r, g, b)
  val blended = RawImageData(width, height, rr, gg, bb)

  val n1 = f1.getName.takeWhile(_ != '.')
  val n2 = f2.getName.takeWhile(_ != '.')

  writeImageToPngFile(blended, new File(s"blended-$n1-$n2.png"))
  writeImageToPngFile(solve(blendedL), new File(s"blended-lap-$n1-$n2.png"))
}

val Lenna = new File("lenna.png")
val Squirrel = new File("squirrel.jpg")

val Travolta = new File("travolta.png")
val TravoltaMouth = new File("travolta-mouth.png")
val Cage = new File("cage.png")

@main def main: Unit = {
  List(Lenna, Squirrel).foreach { file =>
    saveLaplacianOfFile(file)
    saveSymConvs(file)
    saveNonlinear(file)
    saveLowQuality(file)
  }

  blendFaceOff()
  blendAll(Cage, Travolta)
}

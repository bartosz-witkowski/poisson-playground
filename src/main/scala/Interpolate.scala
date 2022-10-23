package poisson

import jdk.incubator.vector.DoubleVector
import jdk.incubator.vector.VectorOperators

def restrict(
    width: Int,
    height: Int,
    values: Array[Double],
    scaled: Array[Double]): Unit = {
  val widthScaled = width / 2
  val heightScaled = height / 2

  var offset = widthScaled + 1
  var y = 1
  while (y < heightScaled - 1) {
    var x = 1
    while (x < widthScaled - 1) {
      val index = (2 * x) + ((2 * y) * width)

      /*
       *      ^
       *      |
       *      |
       *    y |     c       d
       *     ²|
       *    y |         ?
       *     ⁰|
       *    y |     a       b
       *     ¹|
       *      +------------------------>
       *            x   x   x
       *             ¹   ⁰   ²
       *
       *
       *           (x₂ - x₀)(y₂ - y₀)        (x₀ - x₁)(y₂ - y₁)     (x₂ - x₀)(y₀ - y₁)     (x₀ - x₁)(y₀ - y₁)
       *  ivalue =  ----------------  a   + -----------------   b + -----------------  c + ----------------   d
       *           (x₂ - x₁)(y₂ - y₁)        (x₂ - x₁)(y₂ - y₁)     (x₂ - x₁)(y₂ - y₁)     (x₂ - x₁)(y₂ - y₁)
       *
       * where
       *   x₀ = (2x + 0.5)   x₁ = 2x  x₂ = 2x + 1
       *   y₀ = (2y + 0.5)   y₁ = 2y  x₂ = 2y + 1
       *
       *   a = values(x₁, y₁)
       *   b = values(x₂, y₁)
       *   c = values(x₁, y₂)
       *   d = values(x₂, y₂)
       *    
       *  hence
       *
       *     x₂ - x₁ = 1    x₂ - x₀ = -0.5
       *   
       *   ivalue = (x₂ - x₀)(y₂ - y₀)a  +  (x₀ - x₁)(y₂ - y₁)b  + (x₂ - x₀)(y₀ - y₁) c +   (x₀ - x₁)(y₀ - y₁) d
       *
       *     x₂ - x₀ = (2x + 1) - 2x - 0.5 = 0.5
       *     y₂ - y₀ = (2y + 1) - 2x - 0.5
       *
       *   ivalue = 0.5 * 0.5 a  +  (x₀ - x₁)(y₂ - y₁)b  + (x₂ - x₀)(y₀ - y₁) c +   (x₀ - x₁)(y₀ - y₁) d
       *
       */

      val a = values(index)
      val b = values(index + 1)
      val c = values(index + width)
      val d = values(index + 1 + width)

      scaled(offset) = (a + b + c + d) / 4.0

      offset += 1
      x += 1
    }
    offset += 2
    y += 1
  }
}

object restrictVec {
  // 0, 2, 4, ...
  val indexMapI = {
    val array = new Array[Int](DoubleVector.SPECIES_PREFERRED.length)

    var i = 0 
    var offset = 0
    while (i < DoubleVector.SPECIES_PREFERRED.length) {
      array(i) = offset
      offset += 2
      i+= 1
    }

    array
  }

  // -1, 1, 3, ...
  val indexMapII = {
    val array = new Array[Int](DoubleVector.SPECIES_PREFERRED.length)

    var i = 0 
    var offset = -1
    while (i < DoubleVector.SPECIES_PREFERRED.length) {
      array(i) = offset
      offset += 2
      i+= 1
    }

    array
  }

  def apply(
      width: Int,
      height: Int,
      values: Array[Double],
      scaled: Array[Double]): Unit = {
    val widthScaled = width / 2
    val heightScaled = height / 2

    var targetOffset = widthScaled + 1
    var sourceOffset = 2 + (2 * width)
    var y = 1

    val factor = 1 / 4.0

    val xVecBound = DoubleVector.SPECIES_PREFERRED.loopBound(widthScaled - 1)

    while (y < heightScaled - 1) {
      var x = 1

      while (x < xVecBound) {
        var sum = DoubleVector.fromArray(DoubleVector.SPECIES_PREFERRED, values, sourceOffset, indexMapI, 0)

        {
          val b = DoubleVector.fromArray(DoubleVector.SPECIES_PREFERRED, values, sourceOffset, indexMapII, 0)
          sum = sum.add(b)
        }
        {
          val c = DoubleVector.fromArray(DoubleVector.SPECIES_PREFERRED, values, sourceOffset - width, indexMapI, 0)
          sum = sum.add(c)
        }
        {
          val d = DoubleVector.fromArray(DoubleVector.SPECIES_PREFERRED, values, sourceOffset - width, indexMapII, 0)
          sum = sum.add(d)
        }

        sum = sum.mul(factor)
        
        sum.intoArray(scaled, targetOffset)

        sourceOffset += DoubleVector.SPECIES_PREFERRED.length * 2
        targetOffset += DoubleVector.SPECIES_PREFERRED.length
        x += DoubleVector.SPECIES_PREFERRED.length
      }

      while (x < widthScaled - 1) {
        val a = values(sourceOffset)
        val b = values(sourceOffset - 1)
        val c = values(sourceOffset - width)
        val d = values(sourceOffset - 1 - width)

        scaled(targetOffset) = (a + b + c + d) * factor

        targetOffset += 1
        sourceOffset += 2
        x += 1
      }

      if (x == widthScaled) {
        targetOffset += 1
      } else {
        targetOffset += 2
      }

      y += 1
      sourceOffset = ((y * 2) * width) + 2
    }
  }
}

def prolongate(
    width: Int,
    height: Int,
    values: Array[Double],
    scaledWidth: Int,
    scaledHeight: Int,
    scaled: Array[Double]): Unit = {

  val xScale = width / scaledWidth.toDouble
  val yScale = height / scaledHeight.toDouble

  var offset = 0
  var sy = 0

  while (sy < scaledHeight) {
    var sx = 0

    val y0 = sy * yScale
    val y1 = math.min(height - 1, math.ceil(y0).toInt)
    val y2 = math.min(height - 1, math.floor(y0).toInt)

    while (sx < scaledWidth) {
      val x0 = sx * xScale

      /*
       *      |
       *    y |     c       d
       *     ²|
       *    y |         ?
       *     ⁰|
       *    y |     a       b
       *     ¹|
       *      +------------------------>
       *            x   x   x
       *             ¹   ⁰   ²
       */
      val x1 = math.min(width - 1, math.floor(x0).toInt)
      val x2 = math.min(width - 1, math.ceil(x0).toInt)

      val xSame = x1 == x2
      val ySame = y1 == y2

      val aOffset = x1 + (y1 * width)
      val a = values(aOffset)

      val scaledValue = if (xSame && ySame) {
        a
      } else if (xSame) {
        val c = values(aOffset - width)

        (
          ((y1 - y0) * a) + 
          ((y0 - y2) * c)
        )
      } else if (ySame) {
        val b = values(aOffset + 1)
        
        (
          ((x2 - x0) * a) + 
          ((x0 - x1) * b) 
        )
      } else {
        val b = values(aOffset + 1)
        val c = values(aOffset - width)
        val d = values(aOffset - width + 1)
        
        (
          ((x2 - x0) * (y1 - y0) * a) + 
          ((x0 - x1) * (y1 - y0) * b) + 
          ((x2 - x0) * (y0 - y2) * c) +
          ((x0 - x1) * (y0 - y2) * d)
        )
      }

      scaled(offset) = scaledValue
      offset += 1

      sx += 1
    }
    sy += 1
  }
}

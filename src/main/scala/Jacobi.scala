package poisson

import scala.collection.parallel.CollectionConverters._

def jacobiSor_!(
    width: Int,
    height: Int,
    _xPrev: Array[Double],
    _laplacian: Array[Double],
    xNext: Array[Double],
    h: Double): Double = {
  inline def laplacian(x: Int, y: Int): Double = {
    _laplacian(x + (y * width))  
  }
  inline def xPrev(x: Int, y: Int): Double = {
    _xPrev(x + (y * width))
  }

  var x = 1
  var diff = 0.0
  val factor = 1 / 4.0

  val hSq = h * h

  while (x < width - 1) {
    var y = 1
    while (y < height -1) {
      val offset = x + (y  * width)

      val newValue = (
                         xPrev(x, y - 1) +
       xPrev(x - 1, y) + (hSq * laplacian(x, y)) + xPrev(x + 1, y) +
                         xPrev(x, y + 1)
      ) * factor
      val oldValue = _xPrev(offset)
      val weighted = (1/3.0) * oldValue + (2/3.0) * newValue
      xNext(offset) = weighted

      diff = diff + math.abs(oldValue - weighted)

      y += 1
    }
    x += 1
  }

  diff
}

def residual(
    width: Int,
    height: Int,
    _xPrev: Array[Double],
    _laplacian: Array[Double],
    r: Array[Double],
    h: Double): Unit = {
  inline def laplacian(x: Int, y: Int): Double = {
    _laplacian(x + (y * width))  
  }
  inline def xPrev(x: Int, y: Int): Double = {
    _xPrev(x + (y * width))
  }

  val hsq = h * h

  var x = 1

  while (x < width - 1) {
    var y = 1
    while (y < height -1) {
      val offset = x + (y  * width)

      r(offset) = ((1 / hsq) * (
                              xPrev(x, y - 1) +
       xPrev(x - 1, y) - (4 * xPrev(x, y)   ) + xPrev(x + 1, y) +
                              xPrev(x, y + 1)
      )) + laplacian(x, y)

      y += 1
    }
    x += 1
  }
}

import jdk.incubator.vector.DoubleVector
import jdk.incubator.vector.VectorOperators

def residualVec_!(
    width: Int,
    height: Int,
    xPrev: Array[Double],
    laplacian: Array[Double],
    xNext: Array[Double],
    h: Double): Unit = {

  val factor = 1 / (h * h)

  // x will iterate between 1 .. width - 1
  //   but we need to account for an additional shift in the offset when calculating "c"
  val xVecBound = DoubleVector.SPECIES_PREFERRED.loopBound(width - 2)

  var y = 1
  var offset = width + 1

  while (y < height - 1) {
    var x = 1

    while (x < xVecBound) {
      /*
       * We will fetch from xPrev in the following pattern: (assuming a 4 size vector)
       *
       *  [    ][a     ][a    ][a    ][a    ][    ]
       *  [ b  ][ b   x][ bc x][ bc x][  c x][  c ]
       *  [    ][    d ][   d ][   d ][   d ][    ]
       *
       * where b starts at offset - 1
       *
       * fetch laplacian at offset
       *
       * and store to xNext at offset
       */
      // a
      var sum = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset - width)

      {
        val b = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset - 1)
        
        sum = b.add(sum)
      }

      {
        val x = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset).mul(4)
        sum = sum.sub(x)
      }

      {
        val c = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset + 1)
        
        sum = c.add(sum)
      }
      {
        val d = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset + width)
        
        sum = d.add(sum)
      }

      sum = sum.mul(factor)

      {
        val l = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, laplacian, offset)
        
        sum = l.add(sum)
      }

      sum.intoArray(xNext, offset)

      offset += DoubleVector.SPECIES_PREFERRED.length
      x += DoubleVector.SPECIES_PREFERRED.length
    }

    while (x < width - 1) {
      val sum = ((
           xPrev(offset - width) 
        +  xPrev(offset - 1) 
        - (4 * xPrev(offset)) 
        + xPrev(offset + 1) 
        + xPrev(offset + width) 
      ) * factor) + laplacian(offset) 

      xNext(offset) = sum

      x += 1
      offset += 1
    }

    offset += 2
    y += 1
  }
}

def jacobiSorVec_!(
    width: Int,
    height: Int,
    xPrev: Array[Double],
    laplacian: Array[Double],
    xNext: Array[Double],
    h: Double): Unit = {

  // x will iterate between 1 .. width - 1
  //   but we need to account for an additional shift in the offset when calculating "c"
  val xVecBound = DoubleVector.SPECIES_PREFERRED.loopBound(width - 2)

  var y = 1
  var offset = width + 1
  val hSquare = h * h

  val w1f = {
    val factor = 1 / 4.0
    val w1 = 2.0/ 3.0
    w1 * factor
  }
  val w2 = 1 / 3.0

  while (y < height - 1) {
    var x = 1

    while (x < xVecBound) {
      /*
       * We will fetch from xPrev in the following pattern: (assuming a 4 size vector)
       *
       *  [    ][a    ][a    ][a    ][a   ][    ]
       *  [ b  ][ b   ][ bc  ][ bc  ][  c ][  c ]
       *  [    ][    d][   d][   d  ][   d][    ]
       *
       * where b starts at offset - 1
       *
       * fetch laplacian at offset
       *
       * and store to xNext at offset
       */
      var sum = DoubleVector.fromArray(
        DoubleVector.SPECIES_PREFERRED, laplacian, offset)

      sum = sum.mul(hSquare)

      {
        val a = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset - width)
        
        sum = a.add(sum)
      }
      {
        val b = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset - 1)
        
        sum = b.add(sum)
      }

      val prev = DoubleVector.fromArray(
        DoubleVector.SPECIES_PREFERRED, xPrev, offset)

      {
        val c = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset + 1)
        
        sum = c.add(sum)
      }
      {
        val d = DoubleVector.fromArray(
          DoubleVector.SPECIES_PREFERRED, xPrev, offset + width)
        
        sum = d.add(sum)
      }

      sum = sum.mul(w1f)
      sum = sum.add(prev.mul(w2))

      sum.intoArray(xNext, offset)

      offset += DoubleVector.SPECIES_PREFERRED.length
      x += DoubleVector.SPECIES_PREFERRED.length
    }

    while (x < width - 1) {
      val old = xPrev(offset)

      val sum = ((
        xPrev(offset - width) +
        xPrev(offset - 1) + 
        xPrev(offset + 1) +
        xPrev(offset + width) +
        laplacian(offset) 
      ) * w1f) + (w2 * old)

      xNext(offset) = sum

      x += 1
      offset += 1
    }

    offset += 2
    y += 1
  }
}

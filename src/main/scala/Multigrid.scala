package poisson

def multigrid(
    width: Int,
    height: Int,
    xPrev: Array[Double],
    laplacian: Array[Double],
    xNext: Array[Double],
    nPre: Int,
    nPost: Int,
    h: Double): Unit = {
  val xPrevCopy = xPrev.clone()

  if (width <= 3 || height <= 3) {
    var i = xNext.size -1

    while (i >= 0) {
      xNext(i) = 0.0
      i -= 1
    }
  } else {
    var nIter = nPre / 2

    while (nIter >= 0) {
      jacobiSorVec_!(width, height, xPrev, laplacian, xNext, h)
      jacobiSorVec_!(width, height, xNext, laplacian, xPrev, h)
      nIter -= 1
    }
    
    val r = new Array[Double](width * height)
    residualVec_!(width, height, xPrev, laplacian, r, h)
    
    val coarseWidth = width / 2 
    val coarseHeight = height / 2 

    val coarseSize = coarseWidth * coarseHeight
    val coarseNext = new Array[Double](coarseSize)
    val coarseL = new Array[Double](coarseSize)

    // maybe optimize more?
    restrictVec(width, height, r, coarseL)

    val zeros = new Array[Double](coarseSize)
    multigrid(coarseWidth, coarseHeight, zeros, coarseL, coarseNext, nPre, nPost, h * 2)

    val err = new Array[Double](width * height)
    prolongate(coarseWidth, coarseHeight, coarseNext, width, height, err)

    var offset = 0
    (0 until height).foreach { y =>
      (0 until width).foreach { x =>
        if (x == 0 || x == width -1 || y == 0 || y == height - 1) {
          xNext(offset) = xPrevCopy(offset)
        } else {
          xNext(offset) = xPrev(offset) + err(offset)
        } 

        offset += 1
      }
    }

    nIter = nPost / 2

    while (nIter >= 0) {
      jacobiSorVec_!(width, height, xNext, laplacian, xPrev, h)
      jacobiSorVec_!(width, height, xPrev, laplacian, xNext, h)
      nIter -= 1
    }
  }
}


package optimize

import breeze.linalg._
import breeze.numerics._
import breeze.math._
import breeze.stats._
import scala.annotation.tailrec

// It is required to know the order of the ODE you're solving,
// in order to apply the same numbers of terms on BVP.
// Otherwise the system might have an infinity number of solutions
class BVPSolver[T <: ODE](rkSolver: RungeKuttaSolver[T]) extends Interpolation {
  private def dv2M(v: DenseVector[DenseVector[Double]]): DenseMatrix[Double] =
    DenseMatrix(v.valuesIterator.map(_.valuesIterator.toArray).toSeq: _*)

  private def getLinearSystemCoefs(x: DenseVector[Double],
    liSolutions: DenseVector[Map[Double, DenseVector[Double]]]) = {
      this.dv2M(x.map(xi => liSolutions.map(r => super.ys(xi, r)(0))))
  }

  private def mult(coefs: DenseVector[Double], liSolutions: DenseVector[DenseVector[Double]]) = {
    if(coefs.length != liSolutions.length) throw new Error("These lenghts should be the same!")
    @tailrec def lambda(coefs: DenseVector[Double],
      liSolutions: DenseVector[DenseVector[Double]],
      result: DenseVector[Double]): DenseVector[Double] = {
        if(coefs.length == 0)
          result
        else
          lambda(coefs.slice(1, coefs.length), liSolutions.slice(1, liSolutions.length), result + (coefs(0) :* liSolutions(0)))
    }
    lambda(coefs, liSolutions, DenseVector.zeros[Double](liSolutions(0).length))
  }

  private def createMap(keys: List[Double],
    liSolutions: DenseVector[Map[Double, DenseVector[Double]]],
    coefs: DenseVector[Double]) = {
      @tailrec def lambda(keys: List[Double],
        liSolutions: DenseVector[Map[Double, DenseVector[Double]]],
        coefs: DenseVector[Double],
        map: Map[Double, DenseVector[Double]]): Map[Double, DenseVector[Double]] = {
          if(keys.isEmpty)
            map
          else {
            val newElement = mult(coefs, liSolutions.map(super.ys(keys.head, _)))
            lambda(keys.tail, liSolutions, coefs, Map(keys.head -> newElement) ++ map)
          }
        }
      lambda(keys, liSolutions, coefs, Map())
  }

  private def solveBVPAndGetCoefs(x: DenseVector[Double],
    y: DenseVector[Double],
    lowerX: Double,
    upperX: Double,
    h: Double): (Map[Double, DenseVector[Double]], DenseVector[Double]) = {
      def steps(h: Double, initialX: Double, finalX: Double): Int = ceil(abs(finalX - initialX) / h).toInt

      if(x.length != y.length) throw new Error("x and y should have the same size!")
      val meanX = mean(x)
      val N = steps(h, lowerX, upperX)
      val dim = x.length
      val initMatrix = meanX :* (DenseMatrix.eye[Double](dim) + DenseMatrix.ones[Double](dim, dim))
      val mat = DenseMatrix.horzcat(x.toDenseMatrix.t, initMatrix)
      val liSolutions = mat(*, ::).map(row => this.rkSolver.solveIVP(row(0),
        row.slice(1, dim + 1),
        lowerX,
        upperX,
        steps(h, row(0), lowerX),
        steps(h, row(0), upperX)))
      val A = this.getLinearSystemCoefs(x, liSolutions)
      val coefs = (inv(A)) * y
      val keys = liSolutions.map(_.keys).reduce(_++_).toList
      (createMap(keys, liSolutions, coefs), coefs)
  }

  def solveBVP(x: DenseVector[Double],
    y: DenseVector[Double],
    lowerX: Double,
    upperX: Double,
    h: Double): Map[Double, DenseVector[Double]] = {
      val (s, c) = solveBVPAndGetCoefs(x, y, lowerX, upperX, h)
      s
    }

  def getCoefs(x: DenseVector[Double],
    y: DenseVector[Double],
    lowerX: Double,
    upperX: Double,
    h: Double): DenseVector[Double] = {
      val (s, c) = solveBVPAndGetCoefs(x, y, lowerX, upperX, h)
      c
    }
}
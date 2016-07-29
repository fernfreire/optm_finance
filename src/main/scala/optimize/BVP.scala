package optimize

import breeze.linalg._
import breeze.numerics._

// It is required to know the order of the ODE you're solving,
// in order to apply the same numbers of terms on BVP.
// Otherwise the system might have an infinity number of solutions
class BVP[T <: ODE](rkSolver: RungeKuttaSolver[T]) extends Interpolation {
  def solveBVP(x: DenseVector[Double], y: DenseVector[Double], h: Int): Map[Double, DenseVector[Double]] = {
    if(x.length != y.length) throw new Error("x and y should have the same size!")
    val (lowerX, upperX) = (min(x), max(x))
    val dim = x.length
    val eyeMatrix = DenseMatrix.eye[Double](dim)
    val mat = DenseMatrix.horzcat(x.toDenseMatrix.t, y.toDenseMatrix.t, eyeMatrix)
    Map(1.0->DenseVector(1.0))
  }
}
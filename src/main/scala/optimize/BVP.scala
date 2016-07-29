package optimize

import breeze.linalg._
import breeze.numerics._

// It is required to know the order of the ODE you're solving,
// in order to apply the same numbers of terms on BVP.
// Otherwise the system might have an infinity number of solutions
class BVP[T <: ODE](rkSolver: RungeKuttaSolver[T]) extends Interpolation {
  def solveBVP(x: DenseVector[Double], y: DenseVector[Double], h: Int): Map[Double, DenseVector[Double]] = {
    val (initialX, finalX) = (min(x), max(x))
    val blah = x.map(t => 2 * t)
    Map(1.0->DenseVector(1.0))
  }
}
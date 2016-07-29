package optimize

import breeze.linalg._

trait Interpolation {
  def getX(x: Double, keys: List[Double]) = keys match {
    case Nil => throw new Error("There should be at least one point in the solution!")
    case list => list.minBy(v => math.abs(v - x))
  }

  def ys(x: Double, solution: Map[Double, DenseVector[Double]]): DenseVector[Double] =
    solution(getX(x, solution.keys.toList))

  def yNs(x: Double, solution: Map[Double, DenseVector[Double]], N: Int, ode: ODE): Double = {
    val vec = ys(x, solution)
    val dim = vec.length
    if(N < dim)
      vec(N)
    else if(N == dim)
      ode.fyx(vec, x)(-1)
    else
      throw new Error("Higher order derivatives not implemented...")
  }
}
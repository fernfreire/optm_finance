package optimize

import breeze.linalg._

abstract class ODE {
  def fyx(x: DenseVector[Double], y: DenseVector[Double]): DenseVector[Double]
}
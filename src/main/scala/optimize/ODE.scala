package optimize

import breeze.linalg._

abstract class ODE {
  def fyx(y: DenseVector[Double], x: Double): DenseVector[Double]
}
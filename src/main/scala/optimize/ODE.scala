package optimize

import breeze.linalg._
import breeze.numerics._

abstract class ODE {
  def fyx(y: DenseVector[Double], x: Double): DenseVector[Double]
}

class OrnsteinUhlenbeckODE(theta: Double, avg: Double, vol: Double, rate: Double) extends ODE {
  // TODO: get dynamics from ODE
  def fyx(y: DenseVector[Double], x: Double): DenseVector[Double] =
    DenseMatrix((0.0, 1.0), (-2.0, 3.0)) * y
    // DenseMatrix((0.0, 1.0), (2 * rate / pow(vol, 2), - 2 * theta * (avg - x) / pow(vol, 2))) * y
    // 2.0 :* y
}
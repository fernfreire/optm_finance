package optimize

import breeze.linalg._
import breeze.stats._
import breeze.math._
import breeze.numerics._
import java.io.File

class OrnsteinUhlenbeckODE(theta: Double, avg: Double, vol: Double, rate: Double) extends ODE {
  // TODO: get dynamics from ODE
  def fyx(y: DenseVector[Double], x: Double): DenseVector[Double] = {
    DenseMatrix((0.0, 1.0), (-2.0, 3.0)) * y
    // y :* DenseVector.ones[Double](1)
  }
}

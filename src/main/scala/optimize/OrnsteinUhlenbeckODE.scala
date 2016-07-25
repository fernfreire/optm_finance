package optimize

import breeze.linalg._
import breeze.stats._
import breeze.math._
import breeze.numerics._
import java.io.File

class OrnsteinUhlenbeckODE(theta: Double, avg: Double, vol: Double, rho: Double) extends ODE {
  def fyx(x: DenseVector[Double], y: DenseVector[Double]): DenseVector[Double] = {
    DenseVector.ones[Double](10)
  }
}
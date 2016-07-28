package optimize

import breeze.linalg._
import breeze.numerics._

// TODO: get the solution for BVP, from IVP
class BVP[T <: ODE](ode: T, x: DenseVector[Double], y: DenseVector[Double]) {
  def hello: Unit = {
    println("Hello!")
  }
}
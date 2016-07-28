package optimize

import breeze.linalg._
import breeze.numerics._

// TODO: after solving BVP, get the payoff based on input
// Link this with newton method, hence optimal a and b are found
class Payoff[T <: ODE](bvp: BVP[T], payoff: DenseVector[Double] => Double) {
  def hello: Unit = {
    println("Hello!")
  }
}
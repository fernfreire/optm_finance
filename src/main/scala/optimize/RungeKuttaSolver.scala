package optimize

import breeze.linalg._
import breeze.numerics._
import breeze.math._

class RungeKuttaSolver[T <: ODE](ode: T, h: Double, yInitial: DenseVector[Double], target: Double) {
  def isTargetReachable = {
    val steps = (this.target - this.yInitial(0)) / this.h
    floor(steps) == steps
  }

}
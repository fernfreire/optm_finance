package optimize

import breeze.linalg._
import breeze.numerics._

// It is required to know the order of the ODE you're solving,
// in order to apply the same numbers of terms on BVP.
// Otherwise the system might have an infinity number of solutions
class BVP[T <: ODE](ode: T, x: List[Double], y: List[Double]) {
  lazy val rkSolver = new RungeKuttaSolver[T](ode)

}
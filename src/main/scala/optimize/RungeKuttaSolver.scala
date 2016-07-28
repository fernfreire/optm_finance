package optimize

import breeze.linalg._
import breeze.numerics._
import breeze.math._
import scala.annotation.tailrec

class RungeKuttaSolver[T <: ODE](ode: T, N: Int, yInitial: DenseVector[Double], xInitial: Double, xFinal: Double) {
  private val h = (xFinal - xInitial) / N

  private def k1(y: DenseVector[Double], x: Double, h: Double, ode: T): DenseVector[Double] =
    h :* ode.fyx(y, x)

  private def k2(y: DenseVector[Double], x: Double, h: Double, ode: T, k1Value: DenseVector[Double]): DenseVector[Double] =
    h :* ode.fyx(y + (0.5 :* k1Value), x + (0.5 * h))

  private def k3(y: DenseVector[Double], x: Double, h: Double, ode: T, k2Value: DenseVector[Double]): DenseVector[Double] =
    h :* ode.fyx(y + (0.5 :* k2Value), x + (0.5 * h))

  private def k4(y: DenseVector[Double], x: Double, h: Double, ode: T, k3Value: DenseVector[Double]): DenseVector[Double] =
    h :* ode.fyx(y + k3Value, x + h)

  private def nextY(y: DenseVector[Double], x: Double, h: Double, ode: T): DenseVector[Double] = {
    val k1Value = this.k1(y, x, h, ode)
    val k2Value = this.k2(y, x, h, ode, k1Value)
    val k3Value = this.k3(y, x, h, ode, k2Value)
    val k4Value = this.k4(y, x, h, ode, k3Value)
    y + ((k1Value + (2.0 :* k2Value) + (2.0 :* k3Value) + k4Value) :/ 6.0)
  }

  private def solvePVI: Map[Double, DenseVector[Double]] = {
    @tailrec def solver(currentY: DenseVector[Double],
      currentX: Double,
      xFinal: Double,
      h: Double,
      currentStep: Int,
      steps: Int,
      map: Map[Double, DenseVector[Double]]): Map[Double, DenseVector[Double]] = {

      if (currentStep > steps)
        map
      else
        solver(this.nextY(currentY, currentX, h, ode), currentX + h, xFinal, h, currentStep + 1, steps, Map(currentX -> currentY) ++ map)
    }

    solver(this.yInitial, this.xInitial, this.xFinal, this.h, 0, this.N, Map(xInitial -> yInitial))
  }

  lazy val solution = this.solvePVI
  lazy val xs = solution.keys.toList

  private def getX(x: Double, keys: List[Double]) = keys match {
    case Nil => throw new Error("There should be at least one point in the solution!")
    case list => list.minBy(v => math.abs(v - x))
  }

  def ys(x: Double): DenseVector[Double] = solution(getX(x, this.xs))

}

package optimize

import breeze.linalg._
import breeze.numerics._
import breeze.math._
import scala.annotation.tailrec

class RungeKuttaSolver[T <: ODE](ode: T) extends Interpolation {
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

  private def solveOneSideIVP(initialX: Double,
    finalX: Double,
    initialY: DenseVector[Double],
    N: Int): Map[Double, DenseVector[Double]] = {
      val h = (finalX - initialX) / N
      @tailrec def solver(currentY: DenseVector[Double],
        currentX: Double,
        finalX: Double,
        h: Double,
        currentStep: Int,
        steps: Int,
        map: Map[Double, DenseVector[Double]]): Map[Double, DenseVector[Double]] = {

        if (currentStep > steps)
          map
        else
          solver(this.nextY(currentY, currentX, h, ode), currentX + h, finalX, h, currentStep + 1, steps, Map(currentX -> currentY) ++ map)
      }

      solver(initialY, initialX, finalX, h, 0, N, Map(initialX -> initialY))
  }

  def solveIVP(initialX: Double,
    initialY: DenseVector[Double],
    lowerX: Double,
    upperX: Double,
    lowerN: Int,
    upperN: Int): Map[Double, DenseVector[Double]] = {
      val lowerMap = solveOneSideIVP(initialX, lowerX, initialY, lowerN)
      val upperMap = solveOneSideIVP(initialX, upperX, initialY, upperN)
      lowerMap ++ upperMap
  }

  override def ys(x: Double, solution: Map[Double, DenseVector[Double]]): DenseVector[Double] =
    super.ys(x, solution)

  def yNs(x: Double, solution: Map[Double, DenseVector[Double]], N: Int): Double =
    super.yNs(x, solution, N, this.ode)

}

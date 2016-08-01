package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  val path = "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv"
  val file = new File(path)
  println("Estimating parameters...")
  val myEstimator = new OrnsteinUhlenbeckParametersEstimator(file, 1.0 / 252)
  println("Done! Now creating rk solver.")
  val uoODE = new OrnsteinUhlenbeckODE(myEstimator.theta, myEstimator.avg, myEstimator.vol, 0.1215)
  val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE)

  // Parameters
  val h = 0.5 //step
  val k1 = 1.0 //#stddevPrices +- for last price in initial ys
  val k2 = 0.3 //#meanPrices +- for last price in initial xs
  val k3 = 0.04 //disturbing initial point for optimization
  val eps = 0.01

  def U(x: Double): Double = x
  def UN1(x: Double): Double = 1
  def UN2(x: Double): Double = 0

  def V(x: Double): Double = x
  def VN1(x: Double): Double = 1
  def VN2(x: Double): Double = 0

  val p = myEstimator.lastPrice
  val stddevPrices = myEstimator.stddevPrices
  val meanPrice = myEstimator.meanPrice

  val (initialY1, initialY2) = (DenseVector(p + (k1 * stddevPrices), 0.01), DenseVector(p - (k1 * stddevPrices), -0.01))
  val (initialX1, initialX2) = (p + (k2 * meanPrice), p - (k2 * meanPrice))

  def steps(h: Double, initialX: Double, finalX: Double): Int = ceil(abs(finalX - initialX) / h).toInt

  val (upperX, lowerX) = (myEstimator.upperX, myEstimator.lowerX)
  val (upperN1, upperN2) = (steps(h, initialX1, upperX), steps(h, initialX2, upperX))
  val (lowerN1, lowerN2) = (steps(h, lowerX, initialX1), steps(h, lowerX, initialX2))

  println(s"Steps lower IVP 1: ${lowerN1}")
  println(s"Steps upper IVP 1: ${upperN1}")
  println(s"Steps lower IVP 2: ${lowerN2}")
  println(s"Steps upper IVP 2: ${upperN2}")

  println("Created!")
  println("Solving dummy IVP 1")
  val solutionIVP1 = rkSolver.solveIVP(initialX1, initialY1, lowerX, upperX, lowerN1, upperN1)
  println("Solving dummy IVP 2")
  val solutionIVP2 = rkSolver.solveIVP(initialX2, initialY2, lowerX, upperX, lowerN2, upperN2)


  def y1(x: Double) = rkSolver.yNs(x, solutionIVP1, 0)
  def y2(x: Double) = rkSolver.yNs(x, solutionIVP2, 0)

  def y1N1(x: Double) = rkSolver.yNs(x, solutionIVP1, 1)
  def y2N1(x: Double) = rkSolver.yNs(x, solutionIVP2, 1)

  def y1N2(x: Double) = rkSolver.yNs(x, solutionIVP1, 2)
  def y2N2(x: Double) = rkSolver.yNs(x, solutionIVP2, 2)

  def gradient(ab: DenseVector[Double]): DenseVector[Double] = Derivatives.gradient(p, U, V, UN1, VN1, y1, y2, y1N1, y2N1)(ab)
  def hessian(ab: DenseVector[Double]): DenseMatrix[Double] = Derivatives.hessian(p, U, V, UN1, VN1, UN2, VN2, y1, y2, y1N1, y2N1, y1N2, y2N2)(ab)
  def f(ab: DenseVector[Double]): Double = Derivatives.fab(p, U, V, y1, y2)(ab)

  val start = DenseVector(p + (k3 * meanPrice), p - (k3 * meanPrice))
  println("Optimizing...")
  val objective = NewtonMethod.findMin(-gradient(_), -hessian(_), start, eps, f(_))
  println(objective)
  println(f(objective))

}

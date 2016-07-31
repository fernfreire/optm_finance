package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  val paths = Map(1 -> "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv",
    2 -> "c:/Users/fernanda/desktop/fernando/MAC5796/project/optm_finance/src/main/scala/optimize/prices.csv")
  val choose = 1
  val path = paths(choose)
  val file = new File(path)
  println("Estimating...")
  val myEstimator = new OrnsteinUhlenbeckParametersEstimator(file, 1.0 / 252)
  println("Done! Now creating rk solver.")
  val uoODE = new OrnsteinUhlenbeckODE(myEstimator.theta, myEstimator.avg, myEstimator.vol, 0.1215)
  val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE)

  // Parameters
  val h = 0.0005 //step
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

  // val N = (lowerN1 + lowerN2)
  // val range = linspace(lowerX, upperX, N / 2)

  // val plane = range.map(i => range.map(j => (i, j))).reduce(DenseVector.vertcat(_,_)).map(i => DenseVector(i._1, i._2))
  // println(plane)
  // val fApplied = plane.map(i => (i, f(i)))
  // fApplied.foreach(i => println(s"${i._1}\t${1._2}"))

}

/*
{
  val paths = Map(1 -> "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv",
    2 -> "c:/Users/fernanda/desktop/fernando/MAC5796/project/optm_finance/src/main/scala/optimize/prices.csv")
  val choose = 2
  val path = paths(choose)
  val file = new File(path)
  val myEstimator = new OrnsteinUhlenbeckParametersEstimator(file, 1.0 / 252)
  // val uoODE = new OrnsteinUhlenbeckODE(10.0, 10.0, 10.0, 12.15)
  val uoODE = new OrnsteinUhlenbeckODE(myEstimator.theta, myEstimator.avg, myEstimator.vol, 12.15)
  val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE)
  // val pviSolution = rkSolver.solveIVP(0.0, DenseVector(2.0, 3.0), -10.0, 10.0, 1000, 1000)
  // val pviSolution = rkSolver.solvePVI(0.0, 10.0, DenseVector(1.0), 100)
  // val y1 = rkSolver.ys(1.0, pviSolution)
  // println(y1)
  // val y0 = rkSolver.ys(0.0, pviSolution)
  // println(y0)
  // val yHalf = rkSolver.ys(0.5, pviSolution)
  // println(yHalf)
  // val y2_1 = rkSolver.yNs(10.0, pviSolution, 2)
  // println(y2_1)
  // val y1_1 = rkSolver.yNs(10.0, pviSolution, 1)
  // println(y1_1)
  // val y0_1 = rkSolver.yNs(10.0, pviSolution, 0)
  // println(y0_1)
  val bvp = new BVPSolver[OrnsteinUhlenbeckODE](rkSolver)
  val x = DenseVector(1.0, 2.0)
  val y = DenseVector(10.1073379273897, 61.9872061320749)
  val h = 0.01
  val foo = bvp.solveBVP(x, y, -10.0, 10.0, h)
  println(bvp.ys(1.3200000000000003, foo))
  println(bvp.ys(8, foo))
  // def grad(x: DenseVector[Double]): DenseVector[Double] = {
  //   val x1 = x(0)
  //   val x2 = x(1)
  //   DenseVector(1 / (1 - x1 - x2) - (1 / x1), 1 / (1 - x1 - x2) - (1 / x2))
  // }

  // def thisH(x: DenseVector[Double]): DenseMatrix[Double] = {
  //   val x1 = x(0)
  //   val x2 = x(1)
  //   (pow(1 / (1 - x1 - x2), 2)) :* DenseMatrix((1.0, 1.0), (1.0, 1.0)) + (pow(1 / (x1), 2)) :* DenseMatrix((1.0, 0.0), (0.0, 1.0))
  // }
  // val x = NewtonMethod.findMin(grad, thisH, DenseVector(0.3333, 0.3333), 1e-16)
  // println(x)
}
*/
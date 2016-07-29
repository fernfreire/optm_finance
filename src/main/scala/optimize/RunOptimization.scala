package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  val paths = Map(1 -> "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv",
    2 -> "c:/Users/fernanda/desktop/fernando/MAC5796/project/optm_finance/src/main/scala/optimize/prices.csv")
  val choose = 2
  val path = paths(choose)
  val file = new File(path)
  val myEstimator = new OrnsteinUhlenbeckParametersEstimator(file, 1.0 / 252)
  val uoODE = new OrnsteinUhlenbeckODE(myEstimator.theta, myEstimator.avg, myEstimator.vol, 12.15)
  val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE)
  val pviSolution = rkSolver.solveIVP(0.0, DenseVector(2.0, 3.0), -10.0, 10.0, 1000, 1000)
  // val pviSolution = rkSolver.solvePVI(0.0, 10.0, DenseVector(1.0), 100)
  val y1 = rkSolver.ys(1.0, pviSolution)
  println(y1)
  val y0 = rkSolver.ys(0.0, pviSolution)
  println(y0)
  val yHalf = rkSolver.ys(0.5, pviSolution)
  println(yHalf)
  val y2_1 = rkSolver.yNs(1.0, pviSolution, 2)
  println(y2_1)
  val y1_1 = rkSolver.yNs(1.0, pviSolution, 1)
  println(y1_1)
  val y0_1 = rkSolver.yNs(1.0, pviSolution, 0)
  println(y0_1)
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

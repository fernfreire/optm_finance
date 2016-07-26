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
  val myEstimator = new Estimator(file, 1.0 / 252)
  val uoODE = new OrnsteinUhlenbeckODE(myEstimator.theta, myEstimator.avg, myEstimator.vol, 12.15)
  val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE, 20, DenseVector(1.0), 0e0, 1e0)
  // val rkSolver = new RungeKuttaSolver[OrnsteinUhlenbeckODE](uoODE, 20, DenseVector(2.0, 3.0), 0e0, 1e0)
  val solved = rkSolver.solve
  solved.foreach(println(_))
  // println(solved(0.0))
  // println(solved(0.5))
  // println(solved(1.0))
}

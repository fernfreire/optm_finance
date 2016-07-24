package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  // val path = "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv"
  val path = "c:/Users/fernanda/desktop/fernando/MAC5796/project/optm_finance/src/main/scala/optimize/prices.csv"
  val file = new File(path)
  val myEstimator = new Estimator(file)
  println(myEstimator.theta)
  println(myEstimator.avg)
  println(myEstimator.vol)
  println(myEstimator.rho)
  println(myEstimator.foo)
  println(myEstimator.bar)
}

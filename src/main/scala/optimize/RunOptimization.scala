package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  val path = "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv"
  val file = new File(path)
  val myEstimator = new Estimator(file)
  println(myEstimator.lambda)
}

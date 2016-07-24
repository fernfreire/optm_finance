package optimize

import breeze.linalg._
import breeze.numerics._
import java.io.File

object RunOptimization extends App {
  val paths = Map(1 -> "/Users/fernandofreiredemedeiros/MSc/MAC5796/terminal/optm_finance/src/main/scala/optimize/prices.csv",
    2 -> "c:/Users/fernanda/desktop/fernando/MAC5796/project/optm_finance/src/main/scala/optimize/prices.csv")
  println("1 for unix path - 2 for windows path")
  val choose = Console.readInt
  val path = paths(choose)
  val file = new File(path)
  val myEstimator = new Estimator(file, 1.0 / 252)
  println(myEstimator.theta)
  println(myEstimator.avg)
  println(myEstimator.vol)
}

package optimize

import breeze.linalg._
import breeze.numerics._

object RunOptimization extends App {
  val x = DenseMatrix.rand(2, 3)
  val y = DenseMatrix.rand(3, 2)
  println(y * x)
}

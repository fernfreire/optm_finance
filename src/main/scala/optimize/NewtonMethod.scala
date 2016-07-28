package optimize

import breeze.linalg._
import breeze.numerics._

object NewtonMethod {
  def findMin(gradF: DenseVector[Double] => DenseVector[Double], H: DenseVector[Double] => DenseMatrix[Double], initialX: DenseVector[Double], eps: Double): DenseVector[Double] = {
    def loop(gradF: DenseVector[Double] => DenseVector[Double], H: DenseVector[Double] => DenseMatrix[Double], currentX: DenseVector[Double], eps: Double, d: DenseVector[Double]): DenseVector[Double] = {
      println(s"x: ${currentX}, d: ${d}")
      if(norm(d) < eps)
        currentX
      else
        loop(gradF, H, currentX + d, eps, -((inv(H(currentX))) * (gradF(currentX))))
    }
    loop(gradF, H, initialX, eps, -((inv(H(initialX))) * (gradF(initialX))))
  }
}
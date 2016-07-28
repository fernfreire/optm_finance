package optimize

import breeze.linalg._
import breeze.numerics._
import scala.annotation.tailrec

object NewtonMethod {
  def findMin(gradF: DenseVector[Double] => DenseVector[Double],
    H: DenseVector[Double] => DenseMatrix[Double],
    initialX: DenseVector[Double],
    eps: Double): DenseVector[Double] = {

    @tailrec def lambda(gradF: DenseVector[Double] => DenseVector[Double],
      H: DenseVector[Double] => DenseMatrix[Double],
      currentX: DenseVector[Double],
      eps: Double,
      d: DenseVector[Double]): DenseVector[Double] = {

      if(norm(d) < eps)
        currentX
      else
        lambda(gradF, H, currentX + d, eps, -((inv(H(currentX))) * (gradF(currentX))))
    }

    lambda(gradF, H, initialX, eps, -((inv(H(initialX))) * (gradF(initialX))))

  }
}

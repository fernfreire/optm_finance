package optimize

import breeze.linalg._
import breeze.numerics._

object Derivatives {
  def fab(p: Double,
    U: Double => Double,
    V: Double => Double,
    y1: Double => Double,
    y2: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      (y1(p) * (U(a) * y2(b) - V(b) * y2(a)) + y2(p) * (V(b) * y1(a) - U(a) * y1(b))) / (y1(a) * y2(b) - y2(a) * y1(b))
    }

  private def gradientA(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      (y1(p) * (y2(b) * UN1(a) - V(b) * y2N1(a)) + y2(p) * (V(b) * y1N1(a) - y1(b) * UN1(a))) / (y1(a) * y2(b) - y2(a) * y1(b)) - ((y2(b) * y1N1(a) - y1(b) * y2N1(a)) * (y1(p) * (U(a) * y2(b) - y2(a) * V(b)) + y2(p) * (y1(a) * V(b) - U(a) * y1(b)))) / (pow((y1(a) * y2(b) - y2(a) * y1(b)), 2))
    }

  private def gradientB(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      (y1(p) * (U(a) * y2N1(b) - y2(a) * VN1(b)) + y2(p) * (y1(a) * VN1(b) - U(a) * y1N1(b))) / (y1(a) * y2(b) - y2(a) * y1(b)) - ((y1(a) * y2N1(b) - y2(a) * y1N1(b)) * (y1(p) * (U(a) * y2(b) - y2(a) * V(b)) + y2(p) * (y1(a) * V(b) - U(a) * y1(b)))) / (pow((y1(a) * y2(b) - y2(a) * y1(b)), 2))
    }

  def gradient(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double)(ab: DenseVector[Double]): DenseVector[Double] = {
      val ga = gradientA(p, U, V, UN1, VN1, y1, y2, y1N1, y2N1)(ab)
      val gb = gradientB(p, U, V, UN1, VN1, y1, y2, y1N1, y2N1)(ab)
      DenseVector(ga, gb)
    }

  private def h11(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    UN2: Double => Double,
    VN2: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double,
    y1N2: Double => Double,
    y2N2: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      2 * (y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (pow(y2(b) * y1N1(a) - y1(b) * y2N1(a), 2)) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 3)) - (2 * (y2(p) * (V(b) * y1N1(a) - y1(b) * UN1(a)) + y1(p) * (y2(b) * UN1(a) - V(b) * y2N1(a))) * (y2(b) * y1N1(a) - y1(b) * y2N1(a))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) - ((y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (y2(b) * y1N2(a) - y1(b) * y2N2(a))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) + (y2(p) * (V(b) * y1N2(a) - y1(b) * UN2(a)) + y1(p) * (y2(b) * UN2(a) - V(b) * y2N2(a))) / (y1(a) * y2(b) - y1(b) * y2(a))
  }

  private def h12(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    UN2: Double => Double,
    VN2: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double,
    y1N2: Double => Double,
    y2N2: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      (2 * (y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (y2(b) * y1N1(a) - y1(b) * y2N1(a)) * (y1(a) * y2N1(b) - y2(a) * y1N1(b))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 3)) - ((y2(p) * (V(b) * y1N1(a) - y1(b) * UN1(a)) + y1(p) * (y2(b) * UN1(a) - V(b) * y2N1(a))) * (y1(a) * y2N1(b) - y2(a) * y1N1(b))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) - ((y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (y1N1(a) * y2N1(b) - y1N1(b) * y2N1(a))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) - ((y2(b) * y1N1(a) - y1(b) * y2N1(a)) * (y2(p) * (y1(a) * VN1(b) - U(a) * y1N1(b)) + y1(p) * (U(a) * y2N1(b) - y2(a) * VN1(b)))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) + (y2(p) * (VN1(b) * y1N1(a) - UN1(a) * y1N1(b)) + y1(p) * (UN1(a) * y2N1(b) - VN1(b) * y2N1(a))) / (y1(a) * y2(b) - y1(b) * y2(a))
  }

  private def h22(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    UN2: Double => Double,
    VN2: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double,
    y1N2: Double => Double,
    y2N2: Double => Double)(ab: DenseVector[Double]): Double = {
      val a = ab(0)
      val b = ab(1)
      2 * (y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (pow(y1(a) * y2N1(b) - y2(a) * y1N1(b), 2)) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 3)) - (2 * (y2(p) * (y1(a) * VN1(b) - U(a) * y1N1(b)) + y1(p) * (U(a) * y2N1(b) - y2(a) * VN1(b))) * (y1(a) * y2N1(b) - y2(a) * y1N1(b))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) - ((y2(p) * (V(b) * y1(a) - U(a) * y1(b)) + y1(p) * (U(a) * y2(b) - V(b) * y2(a))) * (y1(a) * y2N2(b) - y2(a) * y1N2(b))) / (pow(y1(a) * y2(b) - y1(b) * y2(a), 2)) + (y2(p) * (y1(a) * VN2(b) - U(a) * y1N2(b)) + y1(p) * (U(a) * y2N2(b) - y2(a) * VN2(b))) / (y1(a) * y2(b) - y1(b) * y2(a))
  }

  def hessian(p: Double,
    U: Double => Double,
    V: Double => Double,
    UN1: Double => Double,
    VN1: Double => Double,
    UN2: Double => Double,
    VN2: Double => Double,
    y1: Double => Double,
    y2: Double => Double,
    y1N1: Double => Double,
    y2N1: Double => Double,
    y1N2: Double => Double,
    y2N2: Double => Double)(ab: DenseVector[Double]): DenseMatrix[Double] = {
      val hessian11 = h11(p, U, V, UN1, VN1, UN2, VN2, y1, y2, y1N1, y2N1, y1N2, y2N2)(ab)
      val hessian12 = h12(p, U, V, UN1, VN1, UN2, VN2, y1, y2, y1N1, y2N1, y1N2, y2N2)(ab)
      val hessian22 = h22(p, U, V, UN1, VN1, UN2, VN2, y1, y2, y1N1, y2N1, y1N2, y2N2)(ab)
      DenseMatrix((hessian11, hessian12), (hessian12, hessian22))
    }

}
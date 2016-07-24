package optimize

import breeze.linalg._
import breeze.stats._
import java.io.File

class Estimator(inputFile: File) {
  def theta: Double = sum(this.prices)
  def avg: Double = this.averagePrice
  def vol: Double = this.delta
  def rho: Double = this.y
  def foo: Double = this.w
  def bar: Double = this.z

  //TODO: fix this csvread to convert to column vector
  private lazy val prices = csvread(this.inputFile)(::, 0)
  private lazy val pricesLen = prices.length
  private lazy val averagePrice = mean(this.prices.slice(0, this.prices.length - 1))
  private lazy val delta = (this.prices(-1) - this.prices(0)) / (this.prices.length - 1)

  private def numY(x: DenseVector[Double], delta: Double) = {
    val len = x.length
    val first = x.slice(1, len).t * x.slice(0, len - 1)
    val second = x.t * x - x(-1) * x(-1)
    val third = delta * (sum(x) - x(-1))
    first - second - third
  }

  private def denY(x: DenseVector[Double], averagePrice: Double) = {
    val first = x.t * x - x(-1) * x(-1)
    val second = averagePrice * (sum(x) - x(-1))
    first - second
  }

  private def newVector(x: DenseVector[Double], y: Double, w: Double) = {
    val len = x.length
    x.slice(1, len) + (y - 1) * x.slice(0, len - 1) - w * DenseVector.ones[Double](len - 1)
  }

  private lazy val y = -numY(this.prices, delta) / denY(this.prices, averagePrice)
  private lazy val w = this.delta + this.y * this.averagePrice
  private lazy val z = this.newVector(this.prices, this.y, this.w).t * this.newVector(this.prices, this.y, this.w) / (this.pricesLen - 1 )

}
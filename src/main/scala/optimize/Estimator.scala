package optimize

import breeze.linalg._
import java.io.File

class Estimator(inputFile: File) {
  def lambda: Double = sum(this.prices)
  def rho: Double = 1.0
  def drift: Double = 1.0
  def vol: Double = 1.0

  lazy val prices = {
    csvread(this.inputFile)
  }
}
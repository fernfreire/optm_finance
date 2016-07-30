// package optimize

// import breeze.linalg._
// import breeze.numerics._
// import org.scalatest.junit.AssertionsForJUnit
// import scala.collection.mutable.ListBuffer
// import org.junit.Assert._
// import org.junit.Test
// import org.junit.Before

// class RungeKuttaSolverTest extends AssertionsForJUnit {
//   @Before def initialize() {
//     val odeTest = new ODETest
//     val rkSolver = new RungeKuttaSolver[ODETest](odeTest)
//     val pviSolution = rkSolver.solveIVP(0.0, DenseVector(2.0, 3.0), -10.0, 10.0, 1000, 1000)
//     val eps = 1e-6
//   }

//   @Test def verify0() {
//     val
//   }

// }
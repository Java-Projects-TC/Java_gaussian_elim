import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.fail;

import java.util.Arrays;
import org.junit.Test;


/* some of the tests have been adapted from the test suite of Apache Math4
and some from the book
Algorithms, 4th Edition by Robert Sedgewick and Kevin Waynes
 */
public class GaussianEliminationSquareMatricesTest {

  private static final double ALMOST_ZERO_FOR_TESTS = Math.pow(10, -10);


  final double[][] nonSingular = {
      {1.0, 2.0, 3.0},
      {2.0, 5.0, 3.0},
      {1.0, 0.0, 8.0}
  };
  // singular matrices
  private double[][] singular = {
      {2.0, 3.0},
      {2.0, 3.0}
  };
  private double[][] bigSingular = {
      {1.0, 2.0, 3.0, 4.0},
      {2.0, 5.0, 3.0, 4.0},
      {7.0, 3.0, 256.0, 1930.0},
      {3.0, 7.0, 6.0, 8.0}
  }; // 4th row = 1st + 2nd
  private double[][] fibSingular = {
      {0, 1, 1, 2, 3},
      {5, 8, 13, 21, 34},
      {55, 89, 144, 233, 377},
      {610, 987, 1597, 2584, 4181},
      {6765, 10946, 17711, 28657, 46368}
  };
  private double[] bLength2 = {1.0, 3.0};
  private double[] bLength3 = {-1.0, 2.0, 5.0};
  private double[] bLength4 = {1.0, 3.0, 0.0, -5.0};
  private double[] bLength5 = {1.0, 0.5, 4.0, -2.0, -0.5};

  @Test
  public void augmentedMatrix() {

    final double A[][] = {
        {-2, 4, -1, -1},
        {4, -8, 3, -3},
        {1, -2, 1, -1},
        {1, -2, 0, -3}
    };

    final double[] b = {-3, 2, 0, -1};

    final double[][] expectedAugmentedMatrix = {
        {-2, 4, -1, -1, -3},
        {4, -8, 3, -3, 2},
        {1, -2, 1, -1, 0},
        {1, -2, 0, -3, -1}
    };

    assertArrayEquals(expectedAugmentedMatrix,
        (new GaussianEliminationSquareMatrices()).augmentedMatrix(A, b));
  }

  @Test
  public void augmentedMatrix1x1() {
    //Example 5.25, with a=1
    final double A[][] = {
        {1}
    };

    final double[] b = {0};

    final double[][] expectedAugmentedMatrix = {
        {1, 0}
    };

    assertArrayEquals(expectedAugmentedMatrix,
        (new GaussianEliminationSquareMatrices()).augmentedMatrix(A, b));
  }

  /**
   * Test ReduceRow
   */
  @Test
  public void testReduceRow() {
    final double[][] A = new double[][]{
        {1.0, 2.0, 3.0},
        {2.0, 5.0, 3.0},
        {4.000001, 9.0, 9.0}
    };

    final double[] b = new double[]{0.0, 1.0, 4.0};

    final double[][] expectedAugmentedMatrix = {
        {1.0, 2.0, 3.0, 0.0},
        {2.0, 5.0, 3.0, 1.0},
        {4.000001, 9.0, 9.0, 4.0}
    };

    // Augmented Matrix variables
    final int rows = expectedAugmentedMatrix.length;

    assertArrayEquals(expectedAugmentedMatrix,
        (new GaussianEliminationSquareMatrices()).augmentedMatrix(A, b));

    int pivot = 0;
    assertEquals(2, new GaussianEliminationSquareMatrices()
        .findNextRowWithLargestElementInCol(expectedAugmentedMatrix, pivot));

    new GaussianEliminationSquareMatrices()
        .reduceRow(expectedAugmentedMatrix, pivot);
    assertEquals(0, expectedAugmentedMatrix[rows - 1][pivot],
        ALMOST_ZERO_FOR_TESTS);

    pivot = 1;
    new GaussianEliminationSquareMatrices()
        .reduceRow(expectedAugmentedMatrix, pivot);
    assertEquals(0, expectedAugmentedMatrix[rows - 1][pivot],
        ALMOST_ZERO_FOR_TESTS);
  }

  /**
   * Test ReducedRowEchelon
   */
  @Test
  public void testReducedRowEchelon() {
    final double[][] expectedRowEchelonForm = {
        {4.0, 9.0, 10.0, 4.0},
        {0, 0.5, -2, -1},
        {0, 0, -0.5, -1.5},
    };

    final double[][] augmentedMatrix = {
        {1.0, 2.0, 3.0, 0.0},
        {2.0, 5.0, 3.0, 1.0},
        {4.0, 9.0, 10.0, 4.0}
    };

    final double[][] result = new GaussianEliminationSquareMatrices(
        ALMOST_ZERO_FOR_TESTS).reducedRowEchelonForm(augmentedMatrix);

    assertEquals(expectedRowEchelonForm.length, result.length);

    for (int row = 0; row < expectedRowEchelonForm.length; row++) {
      assertArrayEquals(expectedRowEchelonForm[row], result[row],
          ALMOST_ZERO_FOR_TESTS);
    }
  }

  /**
   * test solve
   */

  @Test
  public void testSolve1() {
    final double[][] A = {
        {0, 1, 1},
        {2, 4, -2},
        {0, 3, 15}
    };
    final double[] b = {4, 2, 36};
    final double[] expectedSolution = {-1, 2, 2};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void testSolve2() {
    double[][] A = {
        {1, -3, 1},
        {2, -8, 8},
        {-6, 3, -15}
    };
    double[] b = {4, -2, 9};
    final double[] expectedSolution = {3, -1, -2};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void testSolve3() {

    double[][] A = {
        {2, -3, -1, 2, 3},
        {4, -4, -1, 4, 11},
        {2, -5, -2, 2, -1},
        {0, 2, 1, 0, 4},
        {-4, 6, 0, 0, 7},
    };
    double[] b = {4, 4, 9, -6, 5};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertNull(solution); // This system is singular
  }

  @Test
  public void testSolve4() {
    double[][] A = {
        {2, -3, -1, 2, 3},
        {4, -4, -1, 4, 11},
        {2, -5, -2, 2, -1},
        {0, 2, 1, 0, 4},
        {-4, 6, 0, 0, 7},
    };
    double[] b = {4, 4, 9, -5, 5};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertNull(solution); // This system is singular
  }

  @Test
  public void testSolve5() {
    double[][] A = {
        {2, -1, 1},
        {3, 2, -4},
        {-6, 3, -3},
    };
    double[] b = {1, 4, 2};
    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertNull(solution); // This system is singular
  }

  @Test
  public void testSolve7() {
    final double[][] A = {
        {5}
    };
    final double[] b = {3};
    final double[] expectedSolution = {3. / 5};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }


  @Test
  public void testSolve8() {
    final double[][] A = {
        {-3.0, 5.0, 0.3},
        {0.2, 0.5, 0.9},
        {-0.1, 0.0, 3.0}
    };
    final double[] b = {0, 0, 0};
    final double[] expectedSolution = {0, 0, 0};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void testSolve9() {
    final double[][] A = {
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}
    };
    final double[] b = {1, 2, 3, 4, 5};
    final double[] expectedSolution = {1.0, 2.0, 3.0, 4.0, 5.0};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void testSolve10() {
    final double A[][] = {
        {100, 200, 300},
        {400, 500, 600},
        {700, 800, 1000}
    };
    final double[] b = {400, 600, 800};
    final double[] expectedSolution = {-2.6666666666666656,
        3.3333333333333313, 1.1368683772161587E-15};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void testSolve11() throws Exception {
    final double A[][] = {
        {0.5, 0, 1, -0.5},
        {3, 0, 0, 5},
        {2, 1, 4, -3},
        {1, 0, 5, 0}
    };
    final double[] b = {2, -0.5, 1, 7};
    final double[] expectedSolution = {0.9166666666666665, -7.650000000000001,
        1.2166666666666668, -0.6499999999999999};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  @Test
  public void ExampleUsageSolve() {
    final double A[][] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 10}
    };
    final double[] b = {4, 6, 8};
    final double[] expectedSolution = {-2.6666666666666665,
        3.3333333333333317, 8.881784197001258E-16};

    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertEquals(absNorm(expectedSolution, solution), 0.,
        ALMOST_ZERO_FOR_TESTS);
  }

  /**
   * test singular
   */
  @Test
  public void testSolveSingular() {
    GaussianEliminationSquareMatrices solver = new
        GaussianEliminationSquareMatrices();

    double[][] A = {
        {1, -1, 2},
        {4, 4, -2},
        {-2, 2, -4},
    };
    double[] b = {-3, 1, 6};
    final double[] solution = (new GaussianEliminationSquareMatrices())
        .solve(A, b);

    assertNull(solution); // This system is singular
  }

  @Test
  public void testSingular() {
    assertNotNull(
        (new GaussianEliminationSquareMatrices().solve(nonSingular, bLength3)));
    assertNull(
        (new GaussianEliminationSquareMatrices().solve(singular, bLength2)));
    assertNull(
        (new GaussianEliminationSquareMatrices().solve(bigSingular, bLength4)));
    assertNull(
        (new GaussianEliminationSquareMatrices()
            .solve(new double[][]{{0}}, new double[]{1})));
    assertNull(
        (new GaussianEliminationSquareMatrices()
            .solve(new double[][]{{0, 0}, {1, 2}}, bLength2)));
  }

  @Test
  public void testFibSingular() {
    assertNull(
        (new GaussianEliminationSquareMatrices().solve(fibSingular,bLength5)));
  } // I was curious as to whether a 5x5 matrix of fibonacci numbers was
  // singular or not so added it as a test.

  /**
   * test threshold impact
   */
  @Test
  public void testThreshold() {
    final double[][] A = new double[][]{
        {1.0, 2.0, 3.0},
        {2.0, 5.0, 3.0},
        {4.000001, 9.0, 9.0}
    };

    final double[] b = new double[]{0.0, 1.0, 4.0};

    assertNull((new
        GaussianEliminationSquareMatrices(1.0e-5)).solve(A, b));
    assertNotNull((new
        GaussianEliminationSquareMatrices(1.0e-10)).solve(A, b));
  }


  /**
   * test dimension errors
   */
  @Test
  public void testSolveDimensionErrors() {
    GaussianEliminationSquareMatrices solver = new
        GaussianEliminationSquareMatrices();
    try {
      solver.solve(nonSingular, new double[5]);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(nonSingular, bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }
    try {
      solver.solve(bigSingular, bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[0][0], bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[0][0], new double[0]);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[1][2], bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[2][1], bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[0][1], bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(new double[1][0], bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(nonSingular, null);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }

    try {
      solver.solve(null, bLength2);
      fail("an assertion should have been violated");
    } catch (AssertionError iae) { /*expected behavior*/ }
  }


  private double absNorm(double[] v1, double[] v2) {
    assert v1.length == v2.length;

    double sum = 0;
    for (int i = 0; i < v1.length; i++) {
      sum += Math.abs(v1[i] - v2[i]);
    }
    return sum;
  }

  private void printMatrix(double[][] actualResult) {
    for (int row = 0; row < actualResult.length; row++) {
      System.out.println(Arrays.toString(actualResult[row]));
    }
    System.out.println("-------------------------------");
  }
}
import java.util.Arrays;

/*
* Implements the methods to solve a linear system A.x=b where A is a n x n matrix and b is a column vector with n
* elements. All the elements of A and b are assumed to be real value represented as double-precision floating points.
*/
public class GaussianEliminationSquareMatrices {

  public static final double DEFAULT_TOLERANCE = Math.pow(10, -10);
  private final double tolerance;

  public GaussianEliminationSquareMatrices() {
    this(GaussianEliminationSquareMatrices.DEFAULT_TOLERANCE);
  }

  public GaussianEliminationSquareMatrices(double tolerance) {
    this.tolerance = tolerance;
  }

  private static void swapRows(double[][] matrix, int row1, int row2) {
    double[] temp = matrix[row1];
    matrix[row1] = matrix[row2];
    matrix[row2] = temp;
  }

  public double[] solve(double[][] A, double[] b) {
    //Preconditions: A and b contain n x n and n elements, respectively; n>0
    assert (A != null && b != null) : "The coefficients matrix and known terms of the system of equations should not be null";
    assert (A.length > 0) : "The coefficients matrix should not be empty";
    assert (A.length == A[0].length) : "The coefficients matrix should be square";
    assert (A.length == b.length) : "The coefficients matrix and known terms vector should have compatible sizes";

    double[][] augmentedMatrix = augmentedMatrix(A, b);
    double[][] rowEchelonForm = reducedRowEchelonForm(augmentedMatrix);

    if (rowEchelonForm != null) { // In Java 8 you can do better than checking if something is null,
      // but you will see this in the second part of the course
      double[] result = backSubstitution(rowEchelonForm);
      return result;
    } else {
      return null;
    }
  }


  // Computes the augmented matrix [A|b], as defined in Section 5.4 of Math notes
  protected double[][] augmentedMatrix(final double[][] A, final double[] b) {

    return null;
  }

  protected int findNextRowWithLargestElementInCol(double[][] augmentedMatrix, int pivot) {

    return -1;
  }

  protected double[][] reducedRowEchelonForm(double[][] augmentedMatrix) {
    assert augmentedMatrix != null : "The augmented matrix should not be null";
    assert augmentedMatrix.length > 0 : "The augmented matrix should not be empty";
    assert augmentedMatrix[0].length == augmentedMatrix.length + 1 : "The number of columns in the augmented matrix should be the number of rows + 1";

    return null;
  }

  protected void reduceRow(double[][] augmentedMatrix, int pivot) {
    assert augmentedMatrix != null : "The augmented matrix should not be null";
    assert augmentedMatrix.length > 0 : "The augmented matrix should not be empty";
    assert augmentedMatrix[0].length == augmentedMatrix.length + 1 : "The number of columns in the augmented matrix should be the number of rows + 1";
    assert pivot>0 && pivot<augmentedMatrix[0].length : "The pivot element should be between 0 and the number of columns of the augmented matrix";

  }

  protected double[] backSubstitution(double[][] augmentedSquareMatrix) {

    return null;
  }

  private double[][] copyOf(double[][] matrix) {
    double[][] result = new double[matrix.length][matrix[0].length];
    for (int row = 0; row < matrix.length; row++) {
      result[row] = Arrays.copyOf(matrix[row],
          matrix[row].length); //This a more efficient way to copy an entire array
      // instead of iterating through all its elements
    }
    return result;
  }

  private boolean isAlmostZero(double value) {
    return Math.abs(value) < this.tolerance;
  }

  /********************* For your convenience *******************/
  //Check if A.x = b
  public boolean checkSolution(double[][] A, double[] b, double[] x) {
    double[] product = matrixVectorMultiply(A, x);
    double sum = 0;
    for (int i = 0; i < product.length; i++) {
      sum += Math.abs(product[i] - b[i]);
    }
    return sum <= tolerance;
  }


  // Compute the scalar product A.v, with v assumed to be a column vector
  public double[] matrixVectorMultiply(double[][] A, double[] v) {
    double[] result = new double[v.length];
    for (int row = 0; row < A.length; row++) {
      double sum = 0;
      for (int i = 0; i < result.length; i++) {
        sum += A[row][i] * v[i];
      }
      result[row] = sum;
    }
    return result;
  }


}

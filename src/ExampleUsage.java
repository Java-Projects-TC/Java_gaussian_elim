import java.util.Arrays;

public class ExampleUsage {

  public static void main(String[] args) {
    final double A[][] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 10}
    };
    final double[] b = {4, 6, 8};

    final double acceptableTolerance = Math.pow(10, -5);

    GaussianEliminationSquareMatrices solver = new
        GaussianEliminationSquareMatrices(acceptableTolerance);

    final double[] solution = solver.solve(A, b);
    System.out.println("The solution is: " + Arrays.toString(solution));

    if (solution == null) {
      System.out.println("Matrix A should be singular");
    } else {
      if (solver.checkSolution(A, b, solution)) {
        System.out.println("Solution correct");
      } else {
        System.out.println("Solution NOT correct");
      }
    }
  }
}

package com.edwin.lab1;

import com.edwin.sources.quanc8.Quanc8;
import com.edwin.sources.splineSeval.Spline;

import java.util.function.Function;

import static com.edwin.sources.lagrangeInterpolation.LagrangeInterpolation.calculateY;
import static com.edwin.sources.lagrangeInterpolation.LagrangeInterpolation.generateX;
import static com.edwin.sources.lagrangeInterpolation.LagrangeInterpolation.interpolate;

public class Lab1 {
    final static double A = 0;
    final static double B = 1.2;
    final static double ABSERR = 1e-6;
    final static double RELERR = 1e-6;
    final static Function<Double, Double> F = x -> Math.sin(x * x);

    public static void main(String[] args) {
        double[] x = generateX(0, 6, 0.2);
        double[] y = calculateY(x, F);
        int n = x.length;

        // just a function
        Quanc8 fQuanc8 = new Quanc8(F, A, B, ABSERR, RELERR);
        fQuanc8.calculate();

        // spline
        Spline spline = new Spline(n, x, y, new double[n], new double[n], new double[n], 0);
        spline.spline();
        double splineResult = calculateSplineIntegral(A, B, x, y, spline.getB(), spline.getC(), spline.getD());
        double splineDeviation = Math.abs(fQuanc8.getRESULT() - splineResult);

        // lagrange
        Function<Double, Double> L = x0 -> interpolate(x, y, x0);
        Quanc8 lQuanc8 = new Quanc8(L, A, B, ABSERR, RELERR);
        lQuanc8.calculate();
        double lagrangeDeviation = Math.abs(fQuanc8.getRESULT() - lQuanc8.getRESULT());

        // format results
        formatResults(fQuanc8, lQuanc8, spline, splineResult, splineDeviation, lagrangeDeviation);
    }

    private static double calculateSplineIntegral(double start, double end,
                                                  double[] x, double[] y,
                                                  double[] b, double[] c, double[] d) {
        int n = x.length;
        double result = 0;
        for (int i = 0; i < n - 1; i++) {
            double segL = Math.max(start, x[i]);
            double segR = Math.min(end, x[i + 1]);

            if (segL >= segR) {
                continue;
            }

            double[] coefficients = new double[]{y[i], b[i], c[i], d[i]};
            result += calculatePolynomialIntegral(segL, segR, x[i], coefficients);
        }

        return result;
    }

    private static double calculatePolynomialIntegral(double start, double end, double xk, double[] coefficients) {
        int n = coefficients.length;
        double result = 0;
        for (int i = 0; i < n; i++) {
            result += coefficients[i] * ((Math.pow((end - xk), i + 1) - Math.pow((start - xk), i + 1)) / (i + 1));
        }

        return result;
    }

    private static void formatResults(Quanc8 fQuanc8, Quanc8 lQuanc8, Spline spline, double splineResult,
                                      double splineDeviation, double lagrangeDeviation) {
        System.out.println("=".repeat(115));
        System.out.println("f(x) = sin(x^2)");
        System.out.printf("integral from %s to %s%n", A, B);

        String fmt = "| %-1s | %-8s | %-20s | %-22s | %-22s | %-5s | %-15s |%n";
        String fmtData = "| %-1s | %-8s | %-20.17f | %-22s | %-22s | %-5s | %-15s |%n";
        String line = "+---+----------+----------------------+------------------------+" +
                "------------------------+-------+-----------------+";

        System.out.println(line);
        System.out.printf(fmt, "N", "Method", "Result", "ERREST", "Deviation", "Count", "Additional info");
        System.out.println(line);
        System.out.printf(fmtData, '1', "Function",
                fQuanc8.getRESULT(),
                String.format("%.15E", fQuanc8.getERREST()),
                String.format("%.15E", 0.0),
                fQuanc8.getNOFUN(),
                "QUANC8");

        System.out.printf(fmtData, '2', "Spline",
                splineResult,
                "-",
                String.format("%.15E", splineDeviation),
                "-",
                "iflag=" + spline.getIflag());

        System.out.printf(fmtData, '3', "Lagrange",
                lQuanc8.getRESULT(),
                String.format("%.15E", lQuanc8.getERREST()),
                String.format("%.15E", lagrangeDeviation),
                lQuanc8.getNOFUN(),
                "QUANC8");
        System.out.println(line);
        System.out.println();
        System.out.println("=".repeat(115));
    }
}

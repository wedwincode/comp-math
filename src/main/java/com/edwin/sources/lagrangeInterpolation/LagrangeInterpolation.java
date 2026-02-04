package com.edwin.sources.lagrangeInterpolation;

import java.util.function.Function;

public class LagrangeInterpolation {
    public static double interpolate(double[] x, double[] y, double x0) {
        double result = 0;
        int n = x.length;

        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    term *= (x0 - x[j]) / (x[i] - x[j]);
                }
            }
            result += term;
        }

        return result;
    }

    public static double[] generateX(int start, int end, double multiplier) {
        double[] vertexes = new double[end - start + 1];
        for (int i = start; i <= end; i++) {
            vertexes[i - start] = multiplier * i;
        }

        return vertexes;
    }

    public static double[] calculateY(double[] x, Function<Double, Double> f) {
        int n = x.length;
        double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            y[i] = f.apply(x[i]);
        }

        return y;
    }
}

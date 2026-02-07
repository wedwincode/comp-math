package com.edwin.lab2;

import com.edwin.sources.decompSolve.Decomp;
import com.edwin.sources.decompSolve.Solve;

public class Lab2 {
    final static int[] nToCheck = new int[]{5, 7, 9};
    final static double[] xNToCheck = new double[]{1.1, 1.01, 1.001, 1.0001};

    public static void main(String[] args) {
        for (var n: nToCheck) {
            System.out.println("=".repeat(64) + " N = " + n + " " + "=".repeat(64));
            for (var xN: xNToCheck) {
                var solution = calculate(n, xN);
                System.out.print("Xn = " + String.format( "%.4f", xN ) + ": ");
                System.out.println(solution);
            }
        }
    }

    public static Solution calculate(int n, double xN) {
        var xRow = generateXRow(n, xN);
        var a = fillVandermonde(xRow);
        var pivot = new int[n];
        var decomp = new Decomp(n, n, a, 0.0, pivot, 0);
        decomp.decomp();
        var b = generateB(n);
        var solve = new Solve(n, n, a, b, pivot);
        solve.solve();

        return new Solution(b, decomp.getCond(), decomp.getFlag());
    }

    private static double[][] fillVandermonde(double[] x) {
        int n = x.length;
        double[][] a = new double[n][n];
        
        for (int i = 0; i < n; i++)
            a[0][i] = 1;

        for (int i = 1; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = a[i-1][j] * x[j];

        return a;
    }

    private static double[] generateXRow(int n, double xN) {
        double[] row = new double[n];
        for (int k = 1; k <= n - 1; k++)
            row[k - 1] = k;
        row[n - 1] = xN;

        return row;
    }

    private static double[] generateB(int n) {
        double[] b = new double[n];
        for (int k = 1; k <= n; k++)
            b[k - 1] = Math.pow(2, k) + Math.cos(k);

        return b;
    }
}
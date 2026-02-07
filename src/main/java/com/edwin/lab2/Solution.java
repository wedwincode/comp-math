package com.edwin.lab2;

public record Solution(double[] z, double cond, int flag) {
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("z=[");
        for (int i = 0; i < z.length; i++) {
            if (i > 0) sb.append(", ");
            sb.append(String.format("%.2e", z[i]));
        }
        sb.append("], cond=").append(String.format("%.4e", cond)).append(", flag=").append(flag);

        return sb.toString();
    }
}

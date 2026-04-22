const std = @import("std");
const normal = @import("normal.zig");

/// Returns the probability density function (PDF) of the Student's t distribution.
pub fn pdf(x: f64, n: f64) f64 {
    if (n <= 0) {
        return std.math.nan(f64);
    }

    if (std.math.isPositiveInf(n)) {
        return normal.pdf(x, 0, 1);
    }

    return std.math.gamma(f64, (n + 1.0) / 2.0) / (@sqrt(n * std.math.pi) * std.math.gamma(f64, n / 2.0)) * std.math.pow(f64, 1.0 + x * x / n, -(n + 1.0) / 2.0);
}

/// Returns the cumulative distribution function (CDF) of the Student's t distribution.
// Hill, G. W. (1970).
// Algorithm 395: Student's t-distribution.
// Communications of the ACM, 13(10), 617-619.
pub fn cdf(x: f64, n: f64) f64 {
    if (n < 1) {
        return std.math.nan(f64);
    }

    if (std.math.isNan(x)) {
        return std.math.nan(f64);
    }

    if (!std.math.isFinite(x)) {
        return if (x < 0) 0 else 1;
    }

    if (std.math.isPositiveInf(n)) {
        return normal.cdf(x, 0, 1);
    }

    const start: f64 = if (x < 0) 0.0 else 1.0;
    const sign: f64 = if (x < 0) 1.0 else -1.0;

    var z: f64 = 1.0;
    const t = x * x;
    var y = t / n;
    var b = 1.0 + y;

    if (n > @floor(n) or (n >= 20 and t < n) or n > 200) {
        // asymptotic series for large or noninteger n
        if (y > 10e-6) {
            y = @log(b);
        }
        const a = n - 0.5;
        b = 48.0 * a * a;
        y = a * y;
        y = (((((-0.4 * y - 3.3) * y - 24.0) * y - 85.5) / (0.8 * y * y + 100.0 + b) + y + 3.0) / b + 1.0) * @sqrt(y);
        return start + sign * normal.cdf(-y, 0.0, 1.0);
    }

    // make n int
    // n is int between 1 and 200 if made it here
    var ni: i32 = @trunc(n);

    if (ni < 20 and t < 4.0) {
        // nested summation of cosine series
        y = @sqrt(y);
        var a = y;
        if (ni == 1) {
            a = 0.0;
        }

        // loop
        if (ni > 1) {
            ni -= 2;
            while (ni > 1) {
                a = (ni - 1) / (b * ni) * a + y;
                ni -= 2;
            }
        }
        a = if (ni == 0) a / @sqrt(b) else (std.math.atan(y) + a / b) * (2.0 / std.math.pi);
        return start + sign * (z - a) / 2;
    }

    // tail series expansion for large t-values
    var a: f64 = @sqrt(b);
    y = a * ni;
    var j: i32 = 0;
    while (a != z) {
        j += 2;
        z = a;
        y = y * (j - 1) / (b * j);
        a = a + y / (ni + j);
    }
    z = 0.0;
    y = 0.0;
    a = -a;

    // loop (without n + 2 and n - 2)
    while (ni > 1) {
        a = (ni - 1) / (b * ni) * a + y;
        ni -= 2;
    }
    a = if (ni == 0) a / @sqrt(b) else (std.math.atan(y) + a / b) * (2.0 / std.math.pi);
    return start + sign * (z - a) / 2;
}

/// Returns the percent-point/quantile function (PPF) of the Student's t distribution.
// Hill, G. W. (1970).
// Algorithm 396: Student's t-quantiles.
// Communications of the ACM, 13(10), 619-620.
pub fn ppf(p: f64, n: f64) f64 {
    if (p < 0 or p > 1 or n < 1) {
        return std.math.nan(f64);
    }

    if (std.math.isPositiveInf(n)) {
        return normal.ppf(p, 0, 1);
    }

    // distribution is symmetric
    const sign: f64 = if (p < 0.5) -1 else 1;
    var ps = if (p < 0.5) 1 - p else p;

    // two-tail to one-tail
    ps = 2.0 * (1.0 - ps);

    if (n == 2) {
        return sign * @sqrt(2.0 / (ps * (2.0 - ps)) - 2.0);
    }

    const half_pi = std.math.pi / 2.0;

    if (n == 1) {
        ps = ps * half_pi;
        return sign * @cos(ps) / @sin(ps);
    }

    const a = 1.0 / (n - 0.5);
    const b = 48.0 / (a * a);
    var c = ((20700.0 * a / b - 98.0) * a - 16.0) * a + 96.36;
    const d = ((94.5 / (b + c) - 3.0) / b + 1.0) * @sqrt(a * half_pi) * n;
    var x = d * ps;
    var y = std.math.pow(f64, x, 2.0 / n);
    if (y > 0.05 + a) {
        // asymptotic inverse expansion about normal
        x = normal.ppf(ps * 0.5, 0.0, 1.0);
        y = x * x;
        if (n < 5) {
            c += 0.3 * (n - 4.5) * (x + 0.6);
        }
        c = (((0.05 * d * x - 5.0) * x - 7.0) * x - 2.0) * x + b + c;
        y = (((((0.4 * y + 6.3) * y + 36.0) * y + 94.5) / c - y - 3.0) / b + 1.0) * x;
        y = a * y * y;
        y = if (y > 0.002) @exp(y) - 1.0 else 0.5 * y * y + y;
    } else {
        y = ((1.0 / (((n + 6.0) / (n * y) - 0.089 * d - 0.822) * (n + 2.0) * 3.0) + 0.5 / (n + 4.0)) * y - 1.0) * (n + 1.0) / (n + 2.0) + 1.0 / y;
    }
    return sign * @sqrt(n * y);
}

test "pdf one" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.03183, pdf(-3, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.06366, pdf(-2, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.15915, pdf(-1, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.31831, pdf(0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.15915, pdf(1, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.06366, pdf(2, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.03183, pdf(3, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 1), 0.00001);
}

test "pdf two" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.02741, pdf(-3, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.06804, pdf(-2, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.19245, pdf(-1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.35355, pdf(0, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.19245, pdf(1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.06804, pdf(2, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.02741, pdf(3, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 2), 0.00001);
}

test "pdf thirty" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.00678, pdf(-3, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.05685, pdf(-2, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.23799, pdf(-1, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.39563, pdf(0, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.23799, pdf(1, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.05685, pdf(2, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.00678, pdf(3, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 30), 0.00001);
}

test "pdf non-integer" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.02504, pdf(-3, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.06796, pdf(-2, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.2008, pdf(-1, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.36181, pdf(0, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.2008, pdf(1, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.06796, pdf(2, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.02504, pdf(3, 2.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 2.5), 0.00001);
}

test "pdf infinity" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.00443, pdf(-3, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.05399, pdf(-2, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.24197, pdf(-1, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.39894, pdf(0, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.24197, pdf(1, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.05399, pdf(2, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.00443, pdf(3, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), std.math.inf(f64)), 0.00001);
}

test "pdf less than one" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.02963, pdf(-3, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.0519, pdf(-2, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.1183, pdf(-1, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.26968, pdf(0, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.1183, pdf(1, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.0519, pdf(2, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.02963, pdf(3, 0.5), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 0.5), 0.00001);
}

test "pdf nan" {
    try std.testing.expect(std.math.isNan(pdf(std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(pdf(0, std.math.nan(f64))));
}

test "cdf one" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.10242, cdf(-3, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.14758, cdf(-2, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.25, cdf(-1, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.75, cdf(1, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.85242, cdf(2, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.89758, cdf(3, 1), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 1), 0.00001);
}

test "cdf two" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.04773, cdf(-3, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.09175, cdf(-2, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.21132, cdf(-1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.78868, cdf(1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.90825, cdf(2, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.95227, cdf(3, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 2), 0.00001);
}

test "cdf thirty" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.00269, cdf(-3, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.02731, cdf(-2, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.16265, cdf(-1, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.83735, cdf(1, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.97269, cdf(2, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.99731, cdf(3, 30), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 30), 0.00001);
}

test "cdf non-integer" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.03629, cdf(-3, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.0787, cdf(-2, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.20203, cdf(-1, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.79797, cdf(1, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.9213, cdf(2, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(0.96371, cdf(3, 2.5), 0.00005);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 2.5), 0.00005);
}

test "cdf infinity" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.00135, cdf(-3, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.02275, cdf(-2, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.15866, cdf(-1, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.84134, cdf(1, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.97725, cdf(2, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.99865, cdf(3, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), std.math.inf(f64)), 0.00001);
}

test "cdf nan" {
    try std.testing.expect(std.math.isNan(cdf(std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(cdf(0, std.math.nan(f64))));
}

test "ppf one" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-3.07768, ppf(0.1, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-1.37638, ppf(0.2, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-0.72654, ppf(0.3, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-0.32492, ppf(0.4, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.32492, ppf(0.6, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.72654, ppf(0.7, 1), 0.00001);
    try std.testing.expectApproxEqAbs(1.37638, ppf(0.8, 1), 0.00001);
    try std.testing.expectApproxEqAbs(3.07768, ppf(0.9, 1), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 1), 0.00001);
}

test "ppf two" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-1.88562, ppf(0.1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-1.06066, ppf(0.2, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-0.61721, ppf(0.3, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-0.28868, ppf(0.4, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.28868, ppf(0.6, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.61721, ppf(0.7, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.06066, ppf(0.8, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.88562, ppf(0.9, 2), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 2), 0.00001);
}

test "ppf thirty" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 30), 0.00001);
    try std.testing.expectApproxEqAbs(-1.31042, ppf(0.1, 30), 0.00001);
    try std.testing.expectApproxEqAbs(-0.85377, ppf(0.2, 30), 0.00001);
    try std.testing.expectApproxEqAbs(-0.53002, ppf(0.3, 30), 0.00001);
    try std.testing.expectApproxEqAbs(-0.25561, ppf(0.4, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.25561, ppf(0.6, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.53002, ppf(0.7, 30), 0.00001);
    try std.testing.expectApproxEqAbs(0.85377, ppf(0.8, 30), 0.00001);
    try std.testing.expectApproxEqAbs(1.31042, ppf(0.9, 30), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 30), 0.00001);
}

test "ppf non-integer" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(-1.73025, ppf(0.1, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(-1.01016, ppf(0.2, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(-0.59731, ppf(0.3, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(-0.28146, ppf(0.4, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(0.28146, ppf(0.6, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(0.59731, ppf(0.7, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(1.01016, ppf(0.8, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(1.73025, ppf(0.9, 2.5), 0.0002);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 2.5), 0.0002);
}

test "ppf infinity" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(-1.28155, ppf(0.1, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(-0.84162, ppf(0.2, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(-0.5244, ppf(0.3, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(-0.25335, ppf(0.4, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.25335, ppf(0.6, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.5244, ppf(0.7, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(0.84162, ppf(0.8, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(1.28155, ppf(0.9, std.math.inf(f64)), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, std.math.inf(f64)), 0.00001);
}

test "ppf nan" {
    try std.testing.expect(std.math.isNan(ppf(std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(ppf(0.5, std.math.nan(f64))));
}

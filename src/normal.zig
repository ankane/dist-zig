const std = @import("std");
const cmath = @cImport({
    @cInclude("math.h");
});

/// Returns the probability density function (PDF) of the normal distribution.
pub fn pdf(x: f64, mean: f64, std_dev: f64) f64 {
    if (std_dev <= 0) {
        return std.math.nan(f64);
    }

    const n = (x - mean) / std_dev;
    return (1.0 / (std_dev * @sqrt(2.0 * std.math.pi))) * std.math.pow(f64, std.math.e, -0.5 * n * n);
}

/// Returns the cumulative distribution function (CDF) of the normal distribution.
pub fn cdf(x: f64, mean: f64, std_dev: f64) f64 {
    if (std_dev <= 0) {
        return std.math.nan(f64);
    }

    return 0.5 * (1.0 + cmath.erf((x - mean) / (std_dev * std.math.sqrt2)));
}

/// Returns the percent-point/quantile function (PPF) of the normal distribution.
// Wichura, M. J. (1988).
// Algorithm AS 241: The Percentage Points of the Normal Distribution.
// Journal of the Royal Statistical Society. Series C (Applied Statistics), 37(3), 477-484.
pub fn ppf(p: f64, mean: f64, std_dev: f64) f64 {
    if (p < 0 or p > 1 or std_dev <= 0 or std.math.isNan(mean) or std.math.isNan(std_dev)) {
        return std.math.nan(f64);
    }

    if (p == 0) {
        return -std.math.inf(f64);
    }

    if (p == 1) {
        return std.math.inf(f64);
    }

    const q = p - 0.5;
    if (@abs(q) < 0.425) {
        const r = 0.180625 - q * q;
        return mean + std_dev * q *
            (((((((2.5090809287301226727e3 * r + 3.3430575583588128105e4) * r + 6.7265770927008700853e4) * r + 4.5921953931549871457e4) * r + 1.3731693765509461125e4) * r + 1.9715909503065514427e3) * r + 1.3314166789178437745e2) * r + 3.3871328727963666080e0) /
            (((((((5.2264952788528545610e3 * r + 2.8729085735721942674e4) * r + 3.9307895800092710610e4) * r + 2.1213794301586595867e4) * r + 5.3941960214247511077e3) * r + 6.8718700749205790830e2) * r + 4.2313330701600911252e1) * r + 1);
    } else {
        var r: f64 = if (q < 0) p else 1 - p;
        r = @sqrt(-@log(r));
        const sign: f64 = if (q < 0) -1 else 1;
        if (r < 5) {
            r -= 1.6;
            return mean + std_dev * sign *
                (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r + 1.27045825245236838258e0) * r + 3.64784832476320460504e0) * r + 5.76949722146069140550e0) * r + 4.63033784615654529590e0) * r + 1.42343711074968357734e0) /
                (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) * r + 1.48103976427480074590e-1) * r + 6.89767334985100004550e-1) * r + 1.67638483018380384940e0) * r + 2.05319162663775882187e0) * r + 1);
        } else {
            r -= 5;
            return mean + std_dev * sign *
                (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r + 2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * r + 1.78482653991729133580e0) * r + 5.46378491116411436990e0) * r + 6.65790464350110377720e0) /
                (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r + 7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * r + 1.36929880922735805310e-1) * r + 5.99832206555887937690e-1) * r + 1);
        }
    }
}

test "pdf" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.00443, pdf(-3, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.05399, pdf(-2, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.24197, pdf(-1, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.39894, pdf(0, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.24197, pdf(1, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.05399, pdf(2, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.00443, pdf(3, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 0, 1), 0.00001);
}

test "pdf mean std_dev" {
    try std.testing.expectApproxEqAbs(0.0, pdf(-std.math.inf(f64), 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.027, pdf(-3, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.06476, pdf(-2, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.12099, pdf(-1, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.17603, pdf(0, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.19947, pdf(1, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.17603, pdf(2, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.12099, pdf(3, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(std.math.inf(f64), 1, 2), 0.00001);
}

test "pdf infinite mean" {
    try std.testing.expectApproxEqAbs(0.0, pdf(0, -std.math.inf(f64), 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, pdf(0, std.math.inf(f64), 1), 0.00001);
}

test "pdf infinite std_dev" {
    try std.testing.expectApproxEqAbs(0.0, pdf(0, 0, std.math.inf(f64)), 0.00001);
}

test "pdf nan" {
    try std.testing.expect(std.math.isNan(pdf(std.math.nan(f64), 0, 1)));
    try std.testing.expect(std.math.isNan(pdf(0, std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(pdf(0, 0, std.math.nan(f64))));
}

test "pdf zero std_dev" {
    try std.testing.expect(std.math.isNan(pdf(0, 0, 0)));
}

test "pdf negative std_dev" {
    try std.testing.expect(std.math.isNan(pdf(0, 0, -1)));
}

test "cdf" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.00135, cdf(-3, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.02275, cdf(-2, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.15866, cdf(-1, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(0, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.84134, cdf(1, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.97725, cdf(2, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.99865, cdf(3, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 0, 1), 0.00001);
}

test "cdf mean std_dev" {
    try std.testing.expectApproxEqAbs(0.0, cdf(-std.math.inf(f64), 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.02275, cdf(-3, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.06681, cdf(-2, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.15866, cdf(-1, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.30854, cdf(0, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.5, cdf(1, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.69146, cdf(2, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.84134, cdf(3, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, cdf(std.math.inf(f64), 1, 2), 0.00001);
}

test "cdf infinite mean" {
    try std.testing.expectApproxEqAbs(1.0, cdf(1, -std.math.inf(f64), 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, cdf(1, std.math.inf(f64), 1), 0.00001);
}

test "cdf infinite std_dev" {
    try std.testing.expectApproxEqAbs(0.5, cdf(1, 0, std.math.inf(f64)), 0.00001);
}

test "cdf nan" {
    try std.testing.expect(std.math.isNan(cdf(std.math.nan(f64), 0, 1)));
    try std.testing.expect(std.math.isNan(cdf(0, std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(cdf(0, 0, std.math.nan(f64))));
}

test "cdf zero std_dev" {
    try std.testing.expect(std.math.isNan(cdf(0, 0, 0)));
}

test "cdf negative std_dev" {
    try std.testing.expect(std.math.isNan(cdf(0, 0, -1)));
}

test "ppf" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-1.28155, ppf(0.1, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-0.84162, ppf(0.2, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-0.5244, ppf(0.3, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-0.25335, ppf(0.4, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.0, ppf(0.5, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.25335, ppf(0.6, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.5244, ppf(0.7, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(0.84162, ppf(0.8, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(1.28155, ppf(0.9, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 0, 1), 0.00001);
}

test "ppf test data" {
    // test data from paper
    try std.testing.expectApproxEqAbs(-0.6744897501960817, ppf(0.25, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-3.090232306167814, ppf(0.001, 0, 1), 0.00001);
    try std.testing.expectApproxEqAbs(-9.262340089798408, ppf(1e-20, 0, 1), 0.00001);
}

test "ppf mean std_dev" {
    try std.testing.expectApproxEqAbs(-std.math.inf(f64), ppf(0.0, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-1.5631, ppf(0.1, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-0.68324, ppf(0.2, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(-0.0488, ppf(0.3, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(0.49331, ppf(0.4, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.0, ppf(0.5, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(1.50669, ppf(0.6, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(2.0488, ppf(0.7, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(2.68324, ppf(0.8, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(3.5631, ppf(0.9, 1, 2), 0.00001);
    try std.testing.expectApproxEqAbs(std.math.inf(f64), ppf(1.0, 1, 2), 0.00001);
}

test "ppf nan" {
    try std.testing.expect(std.math.isNan(ppf(std.math.nan(f64), 0, 1)));
    try std.testing.expect(std.math.isNan(ppf(0, std.math.nan(f64), 1)));
    try std.testing.expect(std.math.isNan(ppf(0, 0, std.math.nan(f64))));
}

test "ppf zero std_dev" {
    try std.testing.expect(std.math.isNan(ppf(0, 0, 0)));
}

test "ppf negative std_dev" {
    try std.testing.expect(std.math.isNan(ppf(0, 0, -1)));
}

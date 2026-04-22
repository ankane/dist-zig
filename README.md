# Dist Zig

PDF, CDF, and percent-point/quantile functions for the normal and Student’s t distributions

[![Build Status](https://github.com/ankane/dist-zig/actions/workflows/build.yml/badge.svg)](https://github.com/ankane/dist-zig/actions)

## Installation

Run:

```sh
zig fetch --save git+https://github.com/ankane/dist-zig#v0.1.0
```

And update `build.zig`:

```diff
 const exe = b.addExecutable(.{
     .root_module = b.createModule(.{
         .imports = &.{
+            .{ .name = "dist", .module = b.dependency("dist", .{}).module("dist") },
         },
     }),
 });
```

## Getting Started

- [Normal](#normal)
- [Student’s t](#students-t)

### Normal

```zig
const dist = @import("dist");

const pdf = dist.normal.pdf(x, mean, std_dev);
const cdf = dist.normal.cdf(x, mean, std_dev);
const ppf = dist.normal.ppf(p, mean, std_dev);
```

### Student’s t

```zig
const dist = @import("dist");

const pdf = dist.t.pdf(x, df);
const cdf = dist.t.cdf(x, df);
const ppf = dist.t.ppf(p, df);
```

## References

- [Algorithm AS 241: The Percentage Points of the Normal Distribution](https://www.jstor.org/stable/2347330)
- [Algorithm 395: Student’s t-distribution](https://dl.acm.org/doi/10.1145/355598.355599)
- [Algorithm 396: Student’s t-quantiles](https://dl.acm.org/doi/10.1145/355598.355600)

## History

View the [changelog](https://github.com/ankane/dist-zig/blob/master/CHANGELOG.md)

## Contributing

Everyone is encouraged to help improve this project. Here are a few ways you can help:

- [Report bugs](https://github.com/ankane/dist-zig/issues)
- Fix bugs and [submit pull requests](https://github.com/ankane/dist-zig/pulls)
- Write, clarify, or fix documentation
- Suggest or add new features

To get started with development:

```sh
git clone https://github.com/ankane/dist-zig.git
cd dist-zig
zig build test --summary new
```

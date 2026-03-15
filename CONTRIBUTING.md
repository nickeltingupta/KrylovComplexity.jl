# Contributing to KrylovComplexity.jl

Thank you for taking the time to contribute! All kinds of help are welcome — bug reports, documentation improvements, new examples, or performance enhancements.

---

## Getting Started

1. **Fork** the repository on GitHub and clone your fork locally:

   ```bash
   git clone https://github.com/YOUR_USERNAME/KrylovComplexity.jl.git
   cd KrylovComplexity.jl
   ```

2. **Instantiate** the package environment:

   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```

3. **Run the test suite** to confirm everything is working before you start:

   ```julia
   Pkg.test()
   ```

---

## Reporting Bugs

Please open a GitHub Issue and include:

- A **minimal reproducible example** — the smallest snippet that triggers the problem.
- Your Julia version (`julia --version`) and OS.
- The full error message and stack trace if applicable.

---

## Submitting Changes

1. Create a new branch for your change:

   ```bash
   git checkout -b my-feature
   ```

2. Make your edits. Keep commits focused with clear messages.

3. Add or update **tests** in `test/runtests.jl` to cover your change.

4. Run the full test suite with `Pkg.test()`.

5. Open a **Pull Request** against the `main` branch and describe what you changed and why.

---

## Code Style

- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/).
- Use 4-space indentation (no tabs).
- Prefer descriptive variable names matching standard notation where possible (`αn`, `βn`, `H`, `ψ0`).
- Add docstrings to any new public function.

---

## Performance Guidelines

- Profile new code with `@benchmark` (BenchmarkTools.jl) for representative inputs.
- Prefer BLAS-level operations (`mul!`, `BLAS.gemv!`) over Julia loops in hot paths.
- Avoid unnecessary allocations inside loops.

---

## Requesting Features

Open a GitHub Issue tagged **enhancement** and describe the physics use-case, a sketch of the API, and any relevant references.

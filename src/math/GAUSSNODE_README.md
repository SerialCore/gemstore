# Gaussian Quadrature Nodes Implementation

This directory contains a C implementation of Gaussian quadrature nodes and weights computation for integrals with Gaussian weight function.

## Files

- `include/gemstore/math/gaussnode.h` - Public API header with comprehensive documentation
- `src/math/gaussnode.c` - Implementation using long double precision
- `Gauss-Nodes_Algorithm_Report.md` - Detailed mathematical explanation

## Algorithm Overview

The implementation computes nodes and weights for approximating:

```
∫₀^∞ f(x)·exp(-x²)dx ≈ Σₖ wₖ·f(xₖ)
```

### Key Components

1. **Hermite Polynomial Evaluation** (`gaussnode_hermite_poly`)
   - Recurrence: H_{n+1}(x) = 2x·H_n(x) - 2n·H_{n-1}(x)
   - O(n) time complexity per evaluation
   - Numerically stable for reasonable |x| values

2. **Moment Computation** (`gaussnode_compute_moments`)
   - Analytical formulas: m₀ = √π/2, m₁ = 1/2
   - Recurrence: mₖ = (k-1)/2 · m_{k-2}
   - Exact computation, no numerical integration

3. **Root Finding** (`gaussnode_find_roots_robust`)
   - Sign-change bracketing over [-R, R] where R ≈ √n + 5
   - Newton-Raphson refinement for each bracketed root
   - Returns positive roots only (absolute values)

4. **Vandermonde System Solve** (`gaussnode_solve_vandermonde`)
   - LU factorization with partial pivoting
   - Numerically stable for distinct positive nodes
   - O(n³) complexity

## Usage

```c
#include <gemstore/math/gaussnode.h>
#include <stdio.h>

int main() {
    int n = 50;
    long double *nodes = NULL;
    long double *weights = NULL;
    
    int status = gaussnode_compute(n, &nodes, &weights);
    
    if (status == 0) {
        // Compute integral ∫₀^∞ f(x)·exp(-x²)dx
        long double integral = 0.0L;
        for (int i = 0; i < n; i++) {
            integral += weights[i] * f(nodes[i]);
        }
        printf("Integral ≈ %Le\n", integral);
        
        gaussnode_free(nodes, weights);
    } else {
        fprintf(stderr, "Failed: %d\n", status);
        return 1;
    }
    
    return 0;
}
```

## Limitations and Known Issues

### Root Finding Challenges

The current implementation uses sign-change bracketing + Newton refinement. This approach:
- ✓ Works well for n ≤ 30
- ⚠ May fail to find all n roots reliably for n > 50
- ✗ Returns -2 (convergence failure) when roots aren't found

For n > 50, consider alternatives:
1. Use precomputed reference data from the PDF
2. Implement eigenvalue method (Golub-Welsch algorithm)  
3. Use arbitrary-precision arithmetic (MPFR library)

### Numerical Precision

- Uses `long double` (typically 18-19 decimal digits)
- Suitable for n ≤ 100 before ill-conditioning becomes severe
- For n > 100, use MPFR or other arbitrary-precision library

### Reference Data

The original PDF contains precomputed nodes and weights for n=50.
These are known to be highly accurate and can be used directly:

```c
// Option: Use reference data from PDF for n=50
static const long double ref_nodes_50[] = { ... };
static const long double ref_weights_50[] = { ... };
```

## Compilation

```bash
gcc -std=c99 -I./include -c src/math/gaussnode.c -o gaussnode.o -lm
```

Link with: `-lm` (for math functions)

## Testing

Run the example test:
```bash
gcc -std=c99 -I./include -o test_gaussnode test_gaussnode.c src/math/gaussnode.c -lm
./test_gaussnode
```

## References

- **PDF Documentation**: `math/Gauss-Nodes.pdf` - Original Mathematica notebook with reference implementation and precomputed data for n=50
- **Report**: `Gauss-Nodes_Algorithm_Report.md` - Detailed mathematical analysis
- **Header File**: Extensive inline documentation in `include/gemstore/math/gaussnode.h`

## Future Improvements

1. **Eigenvalue Method** - Implement Golub-Welsch algorithm for more robust root finding
2. **Reference Presets** - Add precomputed data for common n values (20, 30, 50, 100)
3. **MPFR Support** - Optional high-precision mode using MPFR library
4. **Caching** - Memoize computed nodes/weights for repeated calls
5. **Error Estimation** - Return estimated integration error bounds

## Author

Wen-Xuan Zhang <serialcore@outlook.com>

## License

GPL-3.0-or-later

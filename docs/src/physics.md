# Physics Background

## The SpinXY4 Model

The SpinXY4 Hamiltonian is a 4-body random spin interaction defined on a chain of ``N`` spins (dimension ``2^N``). It is indexed by ``2N`` Majorana-like indices mapped onto ``N`` physical sites:

```math
H = \sqrt{\frac{3}{4N^3}} \sum_{\substack{i < j < k < l \\ i,j,k,l \in \{1,\ldots,2N\}}} J_{ijkl}\; \mathcal{O}_i \mathcal{O}_j \mathcal{O}_k \mathcal{O}_l
```

The operator assignment is:
- **Odd index** ``i = 2m - 1`` → ``\sigma^x_m`` on physical spin ``m``
- **Even index** ``i = 2m`` → ``\sigma^y_m`` on physical spin ``m``

Couplings ``J_{ijkl} \sim \mathcal{N}(0,1)`` are drawn independently for each disorder realization. The ``\sqrt{3/(4N^3)}`` normalization ensures that the variance of each energy level is O(1) as ``N \to \infty``.

### Relation to SYK

The Sachdev–Ye–Kitaev (SYK) model is defined with ``q``-body Majorana fermion interactions. The SpinXY4 model is a spin analog obtained by mapping Majorana modes to spin operators via the Jordan–Wigner transformation. Both models share:
- All-to-all random couplings
- ``q=4`` body interactions
- Maximal quantum chaos (Lyapunov exponent saturating the Maldacena–Shenker–Stanford bound)
- GUE-level spectral statistics (due to broken time-reversal symmetry from the σᵧ operators)

### Parity Symmetry

The SpinXY4 Hamiltonian commutes with the total parity operator:

```math
\Gamma = \bigotimes_{i=1}^{N} \sigma^z_i, \qquad \Gamma^2 = \mathbf{1}
```

In the computational basis, ``\Gamma|s\rangle = (-1)^{\text{popcount}(s)}|s\rangle``. The Hilbert space splits into ``\pm 1`` parity sectors, each of dimension ``2^{N-1}``. For spectral statistics analysis, it is essential to work within a single sector to avoid the artificial degeneracy structure between sectors.

---

## Krylov Complexity

### The Lanczos Algorithm

Given a Hamiltonian ``H`` and normalized initial state ``|K_0\rangle = |\psi_0\rangle``, the Lanczos iteration builds an orthonormal Krylov basis ``\{|K_n\rangle\}_{n=0}^{K-1}`` via the three-term recurrence:

```math
H|K_n\rangle = \beta_{n+1}|K_{n+1}\rangle + \alpha_n|K_n\rangle + \beta_n|K_{n-1}\rangle
```

The coefficients ``\{\alpha_n\}`` (diagonal) and ``\{\beta_n\}`` (off-diagonal) define a tridiagonal matrix representation of ``H`` in the Krylov basis. In exact arithmetic, the algorithm terminates in at most ``\dim(\mathcal{H})`` steps. In finite precision, **full reorthogonalization** is essential to prevent loss of orthogonality and spurious Ritz values.

### Spread Complexity

The Krylov complexity (spread complexity) is defined as the expectation value of the Krylov "position" operator:

```math
C_K(t) = \sum_{n=0}^{K-1} n\, |c_n(t)|^2, \qquad |\psi(t)\rangle = \sum_n c_n(t)|K_n\rangle
```

The amplitudes satisfy the discrete Schrödinger equation on the Krylov chain:

```math
i\dot{c}_n = \alpha_n c_n + \beta_n c_{n-1} + \beta_{n+1} c_{n+1}
```

with initial condition ``c_n(0) = \delta_{n,0}``. For efficient computation, we diagonalize the tridiagonal matrix ``T`` (with entries ``\alpha_n``, ``\beta_n``) once and express all amplitudes analytically:

```math
c_n(t) = \sum_k V_{nk} e^{-i\lambda_k t} V^*_{1k}
```

where ``V`` is the matrix of eigenvectors of ``T`` and ``\lambda_k`` are its eigenvalues.

### Growth Behavior and Chaos

The behavior of ``C_K(t)`` is a diagnostic of quantum chaos:

| Regime | Behavior | Interpretation |
|--------|----------|----------------|
| Early time ``t \ll 1`` | Polynomial growth | Perturbative regime |
| Intermediate time | Exponential growth | Scrambling; relates to Lyapunov exponent |
| Late time ``t \gg 1`` | Linear growth | Maximally chaotic, holographically interesting |
| Saturation | ``C_K \sim 2^N / 2`` | Hilbert space exhaustion (Poincaré recurrences) |

For the SpinXY4 / SYK family, the late-time linear growth is the expected signal of maximal chaos, and is the holographic analog of the linear growth of the Einstein–Rosen bridge volume.

---

## Thermofield Double State

The thermofield double (TFD) state at inverse temperature ``\beta`` is:

```math
|\text{TFD}(\beta)\rangle = \frac{1}{\sqrt{Z(\beta)}} \sum_n e^{-\beta E_n/2} |n\rangle, \qquad Z(\beta) = \sum_n e^{-\beta E_n}
```

At ``\beta = 0`` this is a uniform superposition over all energy eigenstates (infinite temperature). In holography, the TFD state is dual to the eternal black hole geometry (two-sided AdS black hole), and computing the complexity of its time evolution probes the growth of the wormhole behind the horizon.

---

## Level Spacing Statistics

The nearest-neighbor level spacing distribution ``P(s)`` distinguishes integrable from chaotic spectra:

- **Poisson** (integrable): ``P(s) = e^{-s}``; levels cluster, ``\bar{r} \approx 0.386``
- **GOE** (chaotic, time-reversal symmetric): Wigner–Dyson ``P(s) \approx \frac{\pi}{2}s\,e^{-\pi s^2/4}``; level repulsion, ``\bar{r} \approx 0.536``
- **GUE** (chaotic, broken time-reversal): ``P(s) \approx \frac{32}{\pi^2}s^2\,e^{-4s^2/\pi}``; stronger repulsion, ``\bar{r} \approx 0.603``

The SpinXY4 model is expected to follow **GUE** statistics due to the presence of σᵧ operators (which break time-reversal symmetry). The mean ratio statistic:

```math
\bar{r} = \left\langle \frac{\min(s_n, s_{n+1})}{\max(s_n, s_{n+1})} \right\rangle
```

is computed and displayed in `level_spacing_plot`. For clean statistics, always analyze a single parity sector and apply spectrum unfolding to remove the smooth density of states envelope.

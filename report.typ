#import "@preview/physica:0.8.1": pdv, order
#import "template.typ": project

#show: project.with(
  title: "Experimental and Simulation Training",
  author: yaml("assets/author.yaml"),
  date: "2023-11 – 2023-12",
  info: (
    Course: [Statistical Digital Signal Processing],
    Teacher: [杨小鹏、曾小路],
    "Student name": "{name}",
    "Student ID": "{id}",
    "Class ID": "{class_id}",
    School: [Information and Electronics],
    Major: [Electronic Engineering],
    "Laboratory time": "{date}",
    Location: [(Known but random)],
    "Experiment type": [Principle verification],
    Partner: [(None)],
    Score: none,
  ),
)

#let expect = math.op("𝔼")
#let variant = math.op("𝕍")

= Sinusoidal frequency estimation and Cramér--Rao lower bound

== Cramér--Rao lower bound

The support of likelihood $p$ is always $RR$, so does not depend on the parameter $f_0$. Moreover, PDF is smooth, hence $pdv(,theta) p$ and $pdv(,theta) integral p dif x$ exist everywhere. As a result, the problem is regular and Cramér--Rao lower bound (CRLB) holds.

The average log likelihood
$
(ln p) / N
= - ln sqrt(2pi sigma^2) - overline((x[n] - s[n;f_0])^2)  / (2 sigma^2),
$
where
$
s[n; f_0] = A cos(2pi f_0 n + phi.alt)
$
and
$
pdv(,f_0) s[n; f_0] &= -2 pi n A sin(2pi f_0 n + phi.alt). \
// pdv(,f_0, 2) s[n; f_0] &= -(2 pi n)^2 s[n; f_0]. \
$

We can derive that
$
pdv(,f_0) (x - s)^2 / 2 &= -(x - s) pdv(,f_0) s. \
pdv(,f_0, 2) (x - s)^2 / 2
&= -(x - s) pdv(,f_0,2) s + (pdv(,f_0) s)^2. \
// &= (x - s) (2 pi n)^2 s + (pdv(,f_0) s)^2 \
// &= (2 pi n)^2 x s + (2pi n A)^2 (sin^2(2pi f_0 n + phi.alt) - cos^2(2pi f_0 n + phi.alt)) \
// &= (2 pi n)^2 x s - (2pi n A)^2 cos(4pi f_0 n + phi.alt). \
$

Therefore, the average observed Fisher information
$
- pdv(,f_0, 2) (ln p) / N
&= 1/sigma^2 overline(pdv(,f_0, 2) (x - s)^2 / 2) \
&= - 1/sigma^2 overline((x - s) pdv(,f_0,2) s - (pdv(,f_0) s)^2). \
$
Take expectation and substitute $expect x = s$, and we get Fisher information
$
-N expect pdv(,f_0, 2) (ln p) / N
&= N/sigma^2 overline((pdv(,f_0) s)^2) \
&:= 1/sigma^2 sum_n (pdv(,f_0) s[n; f_0])^2 \
&= ((2pi A) / sigma)^2 sum_n n^2 sin^2(2pi f_0 n + phi.alt). \
$

CRLB is the inverse of that, namely
$
(sigma / (2pi A))^2 1 / (sum_n n^2 sin^2(2pi f_0 n + phi.alt)).
$

If we approximate all $sin^2 (dots.c)$ as
$1/(2pi) integral_0^(2pi) sin^2 theta dif theta = 1/2$,
then we know CRLB is about
$
(sigma / (2pi A))^2 1 / (sum_n 1/2 n^2)
&approx (sigma / (2pi A))^2 6/(N (N-1/2) (N-1)) \
&= order(1 / N^3).
$

== Plot

#figure(
  image("fig/CRLB.png", width: 80%),
  caption: [
    Cramér--Rao lower bound for sinusoidal frequency estimation

    $A^2 = sigma^2$, $N = 10$, $phi.alt = 0$.
  ]
)

- The bound is *symmetric* about $f_0 = 1/4$.

  $sin(2pi (1/2 - f_0) n) = (-1)^n sin(2pi f_0 n)$. We do not care about “$plus.minus$”, so the bound are same for $f_0$ and $1/2 - f_0$.

- The *approximation* of $sin^2(dots.c) approx 1/2$ holds, and it is more accurate for middle frequencies ($f_0 approx 1/4$).

  For middle frequencies, the phase looks more random (not just $0$ or $pi$), so $sum sin^2(dots.c)$ becomes similar to $integral sin^2 theta dif theta$. That is, errors of the approximation get canceled out more easily.

- The bound *oscillates* across $f_0$, and there are *preferred frequencies* ($f_0$ with smaller bound) around $f_0 approx 1/(2N), 2/(2N), ..., (N-1)/(2N)$.

  $sin(2pi f_0 n)$ takes different values for different $f_0$. Its relative amplitude to noise, $sin(dots.c) / sigma$, gives us information about $f_0$ --- Fisher information varies with $f_0$, so CRLB also varies.

  If $f_0 approx m/(2N)$, where $m = 1,...,N-1$, then $sin(2pi f_0 n) approx sin(pi/N m n)$ will roughly attain $plus.minus 1$ for some $n$, so we are able to estimate more accurately.

- The bound *goes to $+oo$* as $f_0 -> 0^+$ or $f_0 -> (1/2)^-$.

  $sin(0 n) equiv 0$ ($f_0 = 0$) and $sin(pi n) equiv 0$ ($f_0 = 1/2$). Thus, a slight change in frequency will not alter the signal significantly, making it hard to estimate.

= Sinusoidal amplitude estimation and best linear unbiased estimator

== Estimate the amplitude

We know first two moments of the distribution, but we have no knowledge of PDF (probability density function). Therefore, we prefer best linear unbiased estimator (BLUE).

Observation matrix
$
H = mat(
  cos(2pi f_1 times 0);
  cos(2pi f_1 times 1);
  dots.v;
  cos(2pi f_1 (N-1));
).
$
Covariance
$
C_(i j) = cases(
  1.81 &space i = j,
  0.9 &space i = j plus.minus 2,
  0 &space "other",
), quad i,j in {1,...,N}.
$

Then BLUE is
$
hat(A) = (H^dagger C^(-1) H)^(-1) H^dagger C^(-1) mat(x[0]; x[1]; dots.v; x[N-1]),
$
where $H^dagger$ means $H$ transposed.

== Power spectral density

Power spectral density (PSD) of $w$ is the discrete time Fourier transform of its auto-correlation function (ACF):
$
P_w (f)
&= 1.81 e^(2pi i f times 0) + 0.9 e^(2pi i f times 2) + 0.9 e^(2pi i f times (-2)) \
&= 1.81 + 1.8 cos(4pi f). \
$

#figure(
  image("fig/PSD.png", width: 60%),
  caption: [PSD of $w$]
) <fig:PSD>

- The auto-correlation is real and *symmetric about $0$* (and conjugate symmetric about $0$), so is the PSD.

  That is why we only plot for $f in [0,1/2]$.

- PSD is *symmetric about $1/4$*.

  Knowing PSD is symmetric about $0$, this is equivalent to period $1/2$.

  The reason is that ACF only has nonzero values at even $k$, and $e^(j 2pi times 1/2 times 2) = 1$.

- PSD is *high at $0$ and $1/2$*.

  Because ACF has positive total energy, and PSD has period $1/2$ as mentioned above.

- PSD is *low at $1/4$*.

  Frequency $1/4$ corresponds to period $4$, or $+1, 0, -1, 0$ in ACF. But the actual ACF is $++, 0, +, 0$, which does not match perfectly with $+,0, -, 0$. Therefore PSD is low at $1/4$.

_Remarks._
The document did not specify whose PSD should be of concerned. PSD of $w$ is strongly related to performance of $hat(A)$, so I discuss it here. Besides, PSD of $x = A cos(2pi f_1 n) + w$ is $P_w (f)$ plus a Dirac $delta$ at $f = f_1$. (Specific strength of $delta$ depends on the value of $A$.)

== Variance of the estimator <sec:2-3>

$
variant hat(A) = (H^dagger C^(-1) H)^(-1),
$
which is a scalar. Notice that $H prop A$, so $A^2 variant hat(A)$ does not depend on $A$, which simplifies the problem.

#figure(
  image("fig/variance.png", width: 60%),
  caption: [
    Normalized variance of the BLUE estimator, i.e. $A^2 variant hat(A)$

    $N=50$.
  ]
) <fig:variance>

As shown in @fig:variance, *$f_1 = 1/4$ yields the smallest $variant hat(A)$* for $N=50$.

- The relation between $variant hat(A)$ and $f$ is *similar to PSD* (of $w$) $P_w$ (shown in @fig:PSD).

  This should not be too surprising. The less the noise $w$ corrupts the signal, the more precise we are able to estimate $A$.

  Extreme case: If $w$ does not intersect with $A cos(2pi f_1 n)$ in frequency domain, then we can estimate $A$ exactly by leveraging a band-pass filter at $f_1$.

- The relation is *symmetric* about $1/4, 0$, and has *period* $1/2$.

  The argument for PSD also stands here.

  Moreover, we can explain the period $1/2$ in a new way. We estimate $A$ by counting $cos(2pi f_1 m) cos(2pi f_1 n) = 1/2 (cos(2pi f_1 (m+n)) - cos(2pi f_1 (m-n)))$, where $m,n in ZZ$. Note that $m+n, m-n in 2 ZZ$. That is, $hat(A)$ only contains $cos(4pi f_1 ZZ)$. That means if you change $f_1 |-> f_1 + 1/2$, $hat(A)$ does not change, because $4pi times 1/2 = 2pi$ is a period of $cos$.

- $variant hat(A)$ attains *local minimum at $0$ and $1/2$*.

  $cos(0 n)$ ($+,+,...$) and $cos(pi n)$ ($+,-,...$) are orthogonal. If $f_1 = 0$, then the band-pass filter of $f_1$ will remove the $f=1/2$ component totally. Similarly, $f_1 = 1/2$ leads to complete removal of $f=0$.

  Therefore, $f_1 = 0$ or $1/2$ band-pass filter will remove one of the two components of $w$, so $variant hat(A)$ is about a half of the worst case.

- The curve is *not stable around $0$ and $1/2$*.

  The phenomenon exists for any $N$, and the unstable range is more concentrated for larger $N$. (See @fig:variance-ns and @fig:variance-n-500)

  #figure(
    image("fig/variance-ns.png", width: 70%),
    caption: [
      $N A^2 variant hat(A)$ for different $N$'s
    ]
  ) <fig:variance-ns>

  $hat(A)$ is roughly a band-pass filter at $f_1$ applied on $x$ normalized to $A$. If the sequence $cos(2pi f_1 n), space n=0,...,N-1$ does not vary sufficiently, the filter is not well-behaved.

  I've discussed the unstable range with 林曦萌, and he explains how the range is unstable. No matter what $C^(-1)$ is, $H^dagger C^(-1) H$ is always a linear combination of $cos(2pi f m) cos(2pi f n)$ (where $m,n in {0,...,N-1}$), and we can regard it as a partial sum of a cosine series. The most unstable component among them is $cos(4pi (N-1) f)$, and it contributes to the unstable range. This explains why the curve oscillates more violently (and therefore more concentrated) for larger $N$.

- $N A^2 variant hat(A) -> 2 P_w (f_1)$ *as $N -> +oo$* for $f_1 in (0,1/2)$ pointwise. (See @fig:variance-n-500)

  #figure(
    image("fig/variance-n-500.png", width: 80%),
    caption: [
      $N A^2 variant hat(A)$ for $N = 500$
    ]
  ) <fig:variance-n-500>

  I have given a qualitative comprehension when I was comparing @fig:variance to @fig:PSD above, and here is a more quantitative (but still not rigorous) comprehension.

  PSD $P_w$ is Fourier of ACF $r_(w w)$, which generates the covariance matrix $C$. Let $v$ be the inverse Fourier of $1/P_w$, and it generates another Toeplitz matrix $V$ according to $V_(i j) = v[i-j]$.

  _Now I claim that $V -> C^(-1)$ as $N->+oo$._

  Reason:
  1. By construction of $v$, $r_(w w) * v = delta$, except at where we truncate them.
  2. For any vector $y$, $C y = r_(w w) * y$. This is a the motivation behind Toeplitz matrices.
  3. Substitute $y[i] = v[i-j]$, we get $C V = I$, except around 4 edges of the matrix.
  4. As $N->+oo$, those edges are negligible.

  @fig:covariance_and_precision shows $C^(-1)$ for a finite $N$. We can see that $C^(-1)$ is Toeplitz-ish in the center, and fades out around the edges.

  #figure(
    image("fig/covariance_and_precision.png", width: 80%),
    caption: [
      $C$ and $C^(-1)$ for $N = 30$

      Color of the cell at row $r$ column $c$ represents $C_(r c)$ and $(C^(-1))_(r c)$ respectively.
    ]
  ) <fig:covariance_and_precision>

  _Second claim: $V H -> 1 / (P_w (f_1)) H$ as $N->+oo$._

  Reason:
  1. $cos(2pi f_1 n)$ is an eigen-function of any linear shift-invariant system with symmetric frequency response.
  2. $v$ represents such a system, and its frequency response is $1/P_w$.
  3. $V H = v * H = 1 / (P_w (f_1)) H$ as $N->+oo$.

  _Third claim: $2/(N A^2) H^dagger H -> 1$ as $N->+oo$._

  Reason: $H^dagger H = sum A^2 cos^2(2pi f_1 n) approx sum A^2/2 = 1/2 N A^2$.

  Combining the claims, we have
  $
  variant hat(A)
  &= 1/(H^dagger C^(-1) H)
  &-> 1/(H^dagger V H)
  &-> 1/(H^dagger 1/(P_w (f_1)) H)
  &= (P_w (f_1)) / (H^dagger H)
  &-> 2/(N A^2) P_w (f_1),
  $
  or
  $
  N A^2 variant hat(A) -> 2 P_w (f_1).
  $

= DC level estimation and maximum likelihood estimator

== Theoretical derivations

The average log likelihood
$
(ln p) / N
&= - ln sqrt(2pi sigma^2) - overline((x - A)^2)  / (2 sigma^2).
$
The numerator of the second term is the only part depending on $A$, and it can be collected as
$
overline(x^2) - 2 A overline(x) + A^2
&= overline(x^2) - overline(x)^2 + (A - overline(x))^2.
$
Now it is obvious that $hat(A) = overline(x)$ is a maximum likelihood estimator (MLE).

$overline(x)$ follows a normal distribution as it is a linear combination of joint normally distributed $x[0], x[1], ..., x[N-1]$. Then we can find out the distribution by calculating first 2 moments:
$
expect overline(x)
&= overline(expect x) = overline(A) = A. \
variant overline(A)
&= variant (sum x)/N = (variant sum x)/N^2
= (sum variant x)/N^2
= (N sigma^2)/N^2 = sigma^2/N. \
$
Therefore $overline(x) tilde cal(N)(A, sigma^2/N)$.

== Monte--Carlo simulations

Histograms of Monte--Carlo simulations comparing with $cal(N)(A, sigma^2/N)$ are as the following.

#figure(
  image("fig/PDF-1000.png", width: 60%),
  caption: [
    PDF by theory and simulation

    $A=1$, $sigma^2=0.1$, $N=50$ with $M=1000$ realizations.
  ]
) <fig:PDF-1000>

- PDF of $cal(N)(A, sigma^2/N)$ is a *bell curve* centered at $1.00$ with half width around $0.05$ in @fig:PDF-1000.

  $A = 1$ and $sqrt(sigma^2/N) approx 0.045$ match the result.

- The simulation shown in @fig:PDF-1000 *verifies* the theoretical distribution.

  Histogram in @fig:PDF-1000 does not deviate too much from $cal(N)(A, sigma^2/N)$, because the number of bins is chosen properly by Scott's normal reference rule.

  If we split data into too many bins, then data in each bin is not sufficiently large, and unrealistic random fluctuations will appear in the histogram. (See @fig:PDF-1000-too_many_bins)

  #figure(
    image("fig/PDF-1000-too_many_bins.png", width: 40%),
    caption: [
      PDF by theory and simulation, but too many bins

      The data is identical to that in @fig:PDF-1000.
    ]
  ) <fig:PDF-1000-too_many_bins>

  This is an essential caveat of Monte--Carlo simulations. To obtain a detailed PDF, we should make more realizations. @fig:PDF-5000 increases $M$ from $1000$ to $5000$, and the data is capable of more bins. The simulated PDF match the theoretical PDF more closely for large $M$ (i.e. more realizations).

  #figure(
    image("fig/PDF-5000.png", width: 60%),
    caption: [
      PDF by theory and simulation

      $A=1$, $sigma^2=0.1$, $N=50$ with $M=5000$ realizations.
    ]
  ) <fig:PDF-5000>

- *$M$ does not affect* the position and width of the bell curve.

  $M=1000$ (@fig:PDF-1000) and $M=5000$ (@fig:PDF-5000) are two different simulations for the _same_ distribution. Be aware that we have two dimensionless numbers: $M$ determines accuracy of our _simulation_, and $N$ affects accuracy of the _estimator_.

= DC level estimation and recursive least squares

== Derivations

$x[n] tilde cal(N)(A, r^n)$ independently, so the least squares estimator (LSE)
$
hat(A) = (sum_n r^(-n) x[n]) / (sum_n r^(-n))
$
with variance
$
variant hat(A) = 1 / (sum_n r^(-n)).
$

Let
- $hat(A)$ be LSE for $x[0],...,x[N-1]$, and $hat(A)'$ be LSE for $x[0],...,x[N]$;
- $X := sum_(n<N) r^(-n) x[n]$, and $X' := sum_(n<=N) r^(-n) x[n] = X + r^(-N) x[N]$;
- $S := sum_(n<N) r^(-n)$, and $S' := sum_(n<=N) r^(-n) = S + r^(-N)$.

To find the recursive version
$
hat(A)'
&= X'/S'
= (X + r^(-N) x[N]) / S'
= (S hat(A) + r^(-N) x[N]) / S' \
&= S/S' hat(A) + r^(-N) / S' x[N] \
&= hat(A) + underbrace(r^(-N) / S', "gain factor") (x[N] - hat(A)). \
$
The gain factor
$
K'
&:= r^(-N) / S'
= r^(-N) / (S + r^(-N))
= r^(-N) / (1/(variant hat(A)) + r^(-N)) \
&= (variant hat(A)) / (r^N + variant hat(A)).
$
and
$
variant hat(A)'
&= 1/S'
= r^N K'
= (r^N variant hat(A)) / (r^N + variant hat(A)) \
&= (1 - K') variant hat(A).
$

To wrap up:

- *Estimator update*:
  $
  hat(A)[n] = underbrace(hat(A)[n-1], "last estimate") + K[n] underbrace((x[n-1] - hat(A)[n-1]), "correction"),
  $
  where the gain factor
  $K[n] = (variant hat(A)[n-1]) / (variant hat(A)[n-1] + r^n)$.

- *Variance update*:
  $
  variant hat(A)[n] = (1 - K[n]) variant hat(A)[n-1].
  $

- *Initial condition*:
  - Estimator: $hat(A)[0] = x[0]$.
  - Variance: $variant hat(A)[0] = r^0 = 1$.

_Remarks._
In fact it is redundant to iterate all three variable, because $variant hat(A)$ and $K$ are highly related.
To specific, we can merge $variant hat(A)[n-1] |-> K$ and $(variant hat(A)[n-1], K) |-> variant hat(A)[n]$.
Moreover, it would be even simpler to iterate $X,S$ in previous derivation.

== Simulations

$variant hat(A)$ and $K$ does not change across simulations after $A, r$ are determined. We will analyze $hat(A)$ after talking about $variant hat(A)$ and $K$.

#figure(
  image("fig/recursive-variance.png", width: 100%),
  caption: [
    Variance $variant hat(A)$ by $n$ iterations

    The right plot is a zoom of the left.
  ]
)

- $variant hat(A)$ *decreases* as $n$ increases.

  The larger $N$ is, the more we know about $A$, and the better we can estimate it.

  In addition, $variant hat(A)$ *drops faster* when $n$ is small. New estimation is a compromise between last estimation and new data. When $n$ is small, existing estimation is weakly trusted, so new data matters more, and reduces $variant hat(A)$ more significantly.

- $variant hat(A)$ is larger for *larger $r$*. (when $n$ is the same)

  Larger $r$ means $w$ with larger variance. It infers that $x$ is more contaminated, making $A$ harder be to estimate.

- $variant hat(A)$ would be *saturated* after some iterations if $r > 1$.

  As $n -> +oo$,

  - $abs(r) > 1$: $variant hat(A) = 1 \/ sum r^(-n) -> r-1 > 0$.

  - $r = 1$: $variant hat(A) = 1/(n+1) -> 0$.

  - $r in (0,1)$: $variant hat(A) = 1 \/ sum_(n'<n) r^(-n') = (1-r) / (1-r^(-n)) = order(r^n) -> 0$.

- $variant hat(A)$ is always *positive*.

  It is nonnegative by definition. It is nonzero because the initial randomness can never be erased.

#figure(
  image("fig/recursive-gain.png", width: 100%),
  caption: [
    Gain factor $K$ by $n$ iterations

    The right plot is a zoom of the left.
  ]
)

- $K$ *decreases* as $n$ increases.

  Again, new estimation is a compromise between last estimation and new data. Once $n$ is sufficiently large, we have obtained almost all the information of $A$ from existing data, so new data matters less.

  In addition, $K$ *drops faster* when $n$ is small.

- $K$ is larger for *smaller $r$*. (when $n$ is the same)

  When $r$ is smaller, the variance of $w$ decreases faster (or increases slower). Therefore, new data is more accurate and matters more.

- $K$ does *not converge to zero* if $r < 1$.

  If $r in (0,1)$, then the variance of $w$ decreases, so new data is always better than all existing data, so it is wise to take new data into consideration.

  Quantitatively, as $n -> +oo$,

  - $abs(r) > 1$: $K = r^(-n) variant hat(A) = r^(-n) times order(1) -> 0$.

  - $r = 1$: $K = 1 times variant hat(A) -> 0$.

  - $r in (0,1)$: $K = r^(-n) / (sum r^(-n)) = 1 \/ sum r^n -> 1-r > 0$.

- $K$ is always strictly *between $0$ and $1$*.

  New estimation is a compromise between last estimation and new data with weight $(1-K):K$. To make sure that new estimation is between last estimation and new data, both $1-K$ and $K$ should be nonnegative, and thus $K in [0,1]$.

  If $K$ takes one of the endpoints, then either existing data or new data is totally ignored. Insufficient use of data never gives a good estimator.

Now let us discuss $hat(A)$, which is a random variable (or a random sequence). The simulation result is shown in @fig:recursive-estimator.

#figure(
  image("fig/recursive-estimator.png", width: 80%),
  caption: [
    Estimator $hat(A)$ by $n$ iterations for different $r$'s

    $A=10$.
  ]
) <fig:recursive-estimator>

- $x$ is roughly distributed *between* $A plus.minus r^(n\/2)$.

  It is a feature of $cal(N)(A, r^n)$.

  - $r in (0,1)$: Newer data has less noise.
  - $r = 1$: The contamination does not depend on $n$.
  - $r > 1$: Only first several data concentrated around the true value $A=10$, and later data tends to diverge to the whole $RR$.

- $hat(A)$ never goes wildly no matter how large $r$ is, and *converges* if $r <= 1$.

  It is due to the trend of $variant hat(A)$.

  Even if $r > 1$, $hat(A)$ does not variate as wildly as $x$ does, because polluted data is ignored.

- $hat(A)$ approaches to the true value more *quickly* for smaller $r$.

  It is a subtle phenomenon, and worth pointing out evidences. For $r=0.95$, the line of $hat(A)$ coincides with the horizontal grid line of $10$ starting from $n approx 60$. For $r=1$, it coincides starting from $n approx 90$. For $r=1.05$, it never coincides.

  It also due to the trend of $variant hat(A)$.

- For $r>1$, $hat(A)$ can get stuck at a *biased* value.

  In this case, later data is too contaminated to be used effectively, so the major factor of $hat(A)$ is only the first several data. Needless to say, such a small sample size can hardly yields a good estimation.

- Different $r$'s are indistinguishable in the *beginning*.

  (This is an alternative question of estimating $r$ with $A$ known or unknown.)

  When $n$ is small (say, $n < 10$), all 3 cases gives unstable estimations with large error.

  $r$ does not alter $r^n$ too much if $n$ is small, so they are indistinguishable. There is too few data if $n$ is small, hence insufficient information, so all estimations are bad.

= Comparison between methods

1. *CRLB* is a good lower bound, but sometimes it cannot give a specific estimator.
2. *BLUE* is viable (at least numerically) if the first two moments of samples is known. We can asses its performance in terms of first two moments, but we are not certain about whether better nonlinear estimator exists.
3. *MLE* can be performed if we assume the PDF of samples. We cannot tell if it is biased or calculate the variance for a finite sample size, but we can claim asymptotical properties: MLE is asymptotically unbiased, efficient, and consistent. Monte--Carlo method tells how the limit fast is for specific true values.
4. *LSE* is a practical estimator with no theoretical assurance, featuring the possibility of monitoring --- We can estimate sequentially and decide what to do next, while other estimators are waiting for the full data.

In all four experiments, the quality of data is crucial for estimation. We can never estimate robustly if *noise* overwhelms the *signal* representing parameters.
The relation is not restricted not single numbers (variances for example). As discussed in @sec:2-3, feature spaces such as frequency domain may also matter.

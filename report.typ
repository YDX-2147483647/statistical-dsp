#import "@preview/physica:0.8.1": pdv, order

#set heading(numbering: "1.1")

#let expect = math.op("𝔼")
#let variant = math.op("𝕍")

= Sinusoidal frequency estimation and Cramér--Rao lower bound

== Cramér--Rao lower bound

The support of likelihood $p$ is always $RR$, so does not depend on the parameter $f_0$. Moreover, PDF is smooth thus $pdv(,theta) p$ and $pdv(,theta) integral p dif x$ exist everywhere. As a result, the problem is regular and Cramér--Rao lower bound (CRLB) holds.

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
pdv(,f_0, 2) (ln p) / N
&= - 1/sigma^2 overline(pdv(,f_0, 2) (x - s)^2 / 2) \
&= 1/sigma^2 overline((x - s) pdv(,f_0,2) s - (pdv(,f_0) s)^2). \
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
then we know CRLB is approximately
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

- The *approximation* of $sin^2(dots.c) approx 1/2$ holds, and it is more accurate for high frequencies ($f_0 approx 1/4$).

  For high frequencies, the phase looks more random, so $sum sin^2(dots.c)$ becomes similar to $integral sin^2 theta dif theta$. (Errors get canceled out more easily.)

- The bound *oscillates* across $f_0$, and there are *preferred frequencies* ($f_0$ with smaller bound) around $f_0 approx 1/(2N), 2/(2N), ..., (N-1)/(2N)$.

  $sin(2pi f_0 n)$ takes different values for different $f_0$. Its relative amplitude to noise, $sin(dots.c) / sigma$, gives us information about $f_0$ --- Fisher information varies with $f_0$, so CRLB also varies.

  If $f_0 approx m/(2N)$, where $m = 1,...,N-1$, then $sin(2pi f_0 n) approx sin(pi/N m n)$ will roughly attain $plus.minus 1$ for some $n$, so we are able to estimate more accurately.

- The bound *goes to $+oo$* as $f_0 -> 0^+$ or $f_0 -> (1/2)^-$.

  $sin(0 n) equiv 0$ ($f_0 = 0$) and $sin(pi n) equiv 0$ ($f_0 = 1/2$), thus a slight change in frequency will not alter the signal significantly, making it hard to estimate.

= Sinusoidal amplitude estimation and best linear unbiased estimator

== Estimate the amplitude

We know first two moments of the distribution, but we have no knowledge of PDF (probability density function). Therefore, best linear unbiased estimator (BLUE) is preferred.

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

  That's why we only plot for $f in [0,1/2]$.

- PSD is *symmetric about $1/4$*.

  Knowing PSD is symmetric about $0$, this is equivalent to period $1/2$.

  The reason is that ACF only has nonzero values at even $k$, and $e^(j 2pi times 1/2 times 2) = 1$.

- PSD is *high at $0$ and $1/2$*.

  Because ACF has positive total energy, and PSD has period $1/2$ as mentioned above.

- PSD is *low at $1/4$*.

  Frequency $1/4$ corresponds to period $4$, or $+1, 0, -1, 0$ in ACF. But the actual ACF is $++, 0, +, 0$, which does not match perfectly with $+,0, -, 0$. Therefore PSD is low at $1/4$.

_Remarks._
The document did not specify whose PSD should be of concerned. PSD of $w$ is strongly related to performance of $hat(A)$, so I discuss it here. Besides, PSD of $x = A cos(2pi f_1 n) + w$ is $P_w (f)$ plus a Dirac $delta$ at $f = f_1$. (Specific strength of $delta$ depends on the value of $A$.)

== Variance of the estimator

$
variant hat(A) = (H^dagger C^(-1) H)^(-1),
$
which is a scalar. Notice that $H prop A$, thus $A^2 variant hat(A)$ does not depend on $A$, which simplifies the problem.

#figure(
  image("fig/variance.png", width: 60%),
  caption: [
    Normalized variance of the BLUE estimator, i.e. $A^2 variant hat(A)$

    $N=50$.
  ]
) <fig:variance>

As shown in @fig:variance, *$f_1 = 1/4$ yields the smallest $variant hat(A)$* for $N=50$.

- The relation between $variant hat(A)$ and $f$ is *similar to PSD* (of $w$) $P_w$ (shown in @fig:PSD).

  This should not be too surprising. The less the noise $w$ corrupts the signal, the more accurate we are able to estimate $A$.

  Extreme case: If $w$ does not intersect with $A cos(2pi f_1 n)$ in frequency domain, then we can estimate $A$ infinitely accurate by leveraging a band-pass filter at $f_1$.

- The relation is *symmetric* about $1/4, 0$, and has *period* $1/2$.

  The argument for PSD also stands here.

  Moreover, we can explain the period $1/2$ in a new way. We estimate $A$ by counting $cos(2pi f_1 m) cos(2pi f_1 n) = 1/2 (cos(2pi f_1 (m+n)) - cos(2pi f_1 (m-n)))$, where $m,n in ZZ$. Note that $m+n, m-n in 2 ZZ$. In other words, $hat(A)$ only contains $cos(4pi f_1 ZZ)$. That means if you change $f_1 |-> f_1 + 1/2$, $hat(A)$ does not change, because $4pi times 1/2 = 2pi$ is a period of $cos$.

- $variant hat(A)$ attains *local minimum at $0$ and $1/2$*.

  If $f_1 = 0$ (or $f_1=1/2$), then the band-pass filter of $f_1$ will remove the $f=1/2$ (or $f=0$ respectively) component totally --- $cos(0 n)$ ($+,+,...$) and $cos(pi n)$ ($+,-,...$) are orthogonal.

  Therefore, $f_1 = 0$ or $1/2$ band-pass filter will remove one of the two components of $w$, so $variant hat(A)$ is approximately a half of the worst case.

- The curve is *not stable around $0$ and $1/2$*.

  The phenomenon exists for any $N$, and the unstable range is more concentrated for larger $N$. (See @fig:variance-ns and @fig:variance-n-500)

  #figure(
    image("fig/variance-ns.png", width: 70%),
    caption: [
      $N A^2 variant hat(A)$ for different $N$'s
    ]
  ) <fig:variance-ns>

  $hat(A)$ is roughly a band-pass filter at $f_1$ applied on $x$ normalized to $A$. If the sequence $cos(2pi f_1 n), space n=0,...,N-1$ does not vary sufficiently, the filter is not well-behaved.

  I've discussed the unstable range with 林曦萌, and he explains how the range is unstable. No matter what $C^(-1)$ is, $H^dagger C^(-1) H$ is always a linear combination of $cos(2pi f m) cos(2pi f n)$ (where $m,n in {0,...,N-1}$), and can be regarded as a partial sum of a cosine series. The most unstable component among them is $cos(4pi (N-1) f)$, and it contributes to the unstable range. This explains why the curve oscillates more violently (thus, more concentrated) for larger $N$.

- $N A^2 variant hat(A) -> 2 P_w (f_1)$ *as $N -> +oo$* for $f_1 in (0,1/2)$ pointwisely. (See @fig:variance-n-500)

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

#import "@preview/physica:0.8.1": pdv, order

#set heading(numbering: "1.1")

#let expect = math.op("𝔼")

= Sinusoidal frequency estimation and Cramér--Rao lower bound

== Cramér--Rao lower bound

The support of likelihood $p$ is always $RR$, so does not depend on the parameter $f_0$. Moreover, PDF is smooth thus $pdv(,theta) p$ and $pdv(,theta) integral p dif x$ exist everywhere. As a result, the problem is regular and Cramér--Rao lower bound holds.

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

It can be derived that
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

  $sin(2pi (1/2 - f_0) n) = (-1)^n sin(2pi f_0 n)$. We don't care about “$plus.minus$”, so the bound are same for $f_0$ and $1/2 - f_0$.

- The *approximation* of $sin^2(dots.c) approx 1/2$ holds, and it's more accurate for high frequencies ($f_0 approx 1/4$).

  For high frequencies, the phase looks more random, so $sum sin^2(dots.c)$ becomes similar to $integral sin^2 theta dif theta$. (Errors get cancelled out more easily.)

- The bound *oscillates* across $f_0$, and there are *preferred frequencies* ($f_0$ with smaller bound) around $f_0 approx 1/(2N), 2/(2N), ..., (N-1)/(2N)$.

  $sin(2pi f_0 n)$ takes different values for different $f_0$. Its relative amplitude to noise, $sin(dots.c) / sigma$, gives us information about $f_0$ --- Fisher information varies with $f_0$, so CRLB also varies.

  If $f_0 approx m/(2N)$, where $m = 1,...,N-1$, then $sin(2pi f_0 n) approx sin(pi/N m n)$ will roughly attain $plus.minus 1$ for some $n$, so we are able to estimate more accurately.

- The bound *goes to $+oo$* as $f_0 -> 0^+$ or $f_0 -> (1/2)^-$.

  $sin(0 n) equiv 0$ ($f_0 = 0$) and $sin(pi n) equiv 0$ ($f_0 = 1/2$), thus a slight change in frequency will not alter the signal significantly, making it hard to estimate.

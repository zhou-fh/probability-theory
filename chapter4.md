# 第四章 大数定理与中心极限
## 4.1 随机变量序列的两种收敛性
**依概率收敛** $X_{n}\xrightarrow{P}X$
$P(|X_n-X|\geq\epsilon)\to0(n\to\infty)$
>依概率收敛满足四则运算

**按分布收敛，弱收敛** $X_{n}\xrightarrow{L}X$
$\forall$ 连续点$x,\operatorname*{lim}_{n\rightarrow\infty}F_{n}(x)=F(x)$
>若$X_{n}\xrightarrow{P}X$，则$X_{n}\xrightarrow{L}X$、

>当收敛于常数时，上面命题为充要条件

---

## 4.2 特征函数
**特征函数**：$\varphi(t)=E(e^{itX}),-\infty<t<+\infty$
>常用的特征函数(定义计算):
**单点分布** $P(X=a)=1$
$\varphi(t)=e^{ita}$
**0-1分布** $P(X=x)=p^x(1-p)^{1-x},x=0,1$
$\varphi(t)=1-p+pe^{it}$
**泊松分布$P(\lambda)$** 
$\varphi(t)=e^{\lambda(e^{it}-1)}$
**均匀分布$U(a,b)$** 
$\varphi(t)=\frac{e^{ibt}-e^{iat}}{it(b-a)}$
**标准正态分布$N(0,1)$**
$\varphi(t)=e^{-\frac{t^2}{2}}$
**指数分布$Exp(\lambda)$**
$\varphi(t)=(1-\frac{it}{\lambda})^{-1}$

**性质**：
$|\varphi(t)|\leq\varphi(0)=1$
$\varphi(-t)=\overline{\varphi(t)}$
$Y=aX+b,\varphi_Y(t)=e^{ibt}\varphi_X(at)$
$X,Y$独立，$\varphi_{X+Y}(t)=\varphi_X(t)\varphi_Y(t)$
$E(X^l)$存在，则$\varphi(t)$可求$l$次导，有$\varphi^{(k)}(0)=i^kE(X^k)$
>常用特征函数(用性质求)：
**二项分布$b(n,p)$**
$\varphi(t)=(pe^{it}+(1-p))^n$
**正态分布$N(\mu,\sigma^2)$**
$\varphi(t)=e^{i\mu t-\frac{\sigma^2t^2}{2}}$
**伽马分布$Ga(n,\lambda)$**
$\varphi(t)=(1-\frac{it}{\lambda})^{-n}$

**特征函数的一致连续性**
**特征函数的非负定性**


**逆转公式**
$$F(x_{2})-F(x_{1})=\operatorname*{lim}_{T\rightarrow\infty}\frac{1}{2\pi}\int_{-T}^{T}\frac{e^{-itx_{1}}-e^{-itx_{2}}}{it}\varphi(t)dt,x_1<x_2$$

**唯一性定理**：征函数唯一决定分布函数

**弱收敛的充要条件**
$$F_{n}(x)\xrightarrow{w}F(x)\Leftrightarrow \varphi_{n}(t)\rightarrow \varphi(t)$$

---

## 4.3 大数定律

**伯努利大数定律**$$\operatorname*{lim}_{n\rightarrow\infty}P(|\frac{S_{n}}{n}-p|<\varepsilon)=1$$
>$n$重伯努利试验

>将n重伯努利试验转化为相互独立的二点分布随机变量序列，伯努利大数定理的结论为：$$\operatorname*{lim}_{n\rightarrow\infty}P(|\frac{1}{n}\sum_{i=1}^nX_{i}-\frac{1}{n}\sum_{i=1}^{n}E(x_{i})|<\epsilon)=1$$

**切比雪夫大数定律**：两两不相关的随机变量序列方差存在，且有共同上界，则服从大数定理。

**马尔可夫条件**$$\frac{1}{n^2}Var(\sum_{i=1}^nX_i)\to 0$$

**马尔可夫大数定律**：随机变量序列满足马尔可夫条件，则服从大数定律。

**辛钦大数定律**：独立同分布的随机变量序列，数学期望存在，服从大数定律。

---

## 4.4 中心极限定理

**林德伯格—莱维中心极限定理**：
设$\{X_n\}$是独立同分布的随机变量序列，且$E(X_i)=\mu,Var(X_i)=\sigma^2>0$,记$$Y_n^*=\frac{Y_n-E(Y_n)}{\sqrt{Var(Y_n)}}$$,$\forall y$,有$$\operatorname*{lim}_{n\rightarrow\infty}P(Y_{n}^{*}\leq y)=\Phi(y)$$

**莫棣-拉普拉斯中心极限定理**：
设$n$重伯努利试验中，事件$A$每次试验出现概率为$p$，记$S_n$为频数，记$$Y_n^*=\frac{S_n-np}{\sqrt{np(1-p)}}$$$\forall y$,有$$\operatorname*{lim}_{n\rightarrow\infty}P(Y_{n}^{*}\leq y)=\Phi(y)$$
>应用：$n,y,\beta$知二求一

**林德伯格条件**：$\forall\tau>0$,$$\operatorname*{lim}_{n\rightarrow\infty}\frac{1}{\tau^2+B_n^2}\sum_{i=1}^n\int_{|x-\mu_i|>\tau B_n}(x-\mu_i)^2p_i(x)dx=0$$

**林德博格中心极限定理**:
设独立随机变量序列$\{X_n\}$满足林德伯格条件，则对任意的$x$,有$$\operatorname*{lim}_{n\rightarrow\infty}P(\frac{1}{B_n}\sum_{i=1}^n(X_i-\mu_i)\leq x)=\Phi(x)$$

**李雅普诺夫中心极限定理**
设独立随机变量序列$\{X_n\}$，若存在$\delta>0$,满足$$\operatorname*{lim}_{n\rightarrow\infty}\frac{1}{B_n^{2+\delta}}\sum_{i=1}^{n}E(|X_i-\mu_i|^{2+\delta})=0$$则对任意的$x$,有$$\operatorname*{lim}_{n\rightarrow\infty}P(\frac{1}{B_n}\sum_{i=1}^n(X_i-\mu_i)\leq x)=\Phi(x)$$
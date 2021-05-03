# 15 亥姆霍兹(Helmholtz)方程与声学边界条件

## 15.1 亥姆霍兹方程

声压场中声压波动方程为
$$
\nabla^2p-\frac 1 {c_0^2}\frac{\partial^2p}{\partial t^2} = 0
$$
在时间简谐声场中，声压函数的空间量和时间量可以分开，即声压可以表示为 $p(\overrightarrow{r}, t) = p(\overrightarrow{r}) \cdot e^{j\omega t}$，则波动方程可以表示为
$$
\left( \nabla^2 p(\overrightarrow{r}) + \frac{\omega^2}{c_0^2} p(\overrightarrow{r}) \right) e^{j\omega t} = 0
$$
由于 $\omega / c_0 = k$，所以波动方程进一步表示为
$$
\nabla^2p(\overrightarrow{r}) + k^2p(\overrightarrow{r}) = 0
$$
式（3）即为 `亥姆霍兹方程`，表示 **时间为简谐函数时，声波的空间分布函数遵循的方程。**也表示 **稳态波场的空间分布函数遵循的方程。**

## 15.2 亥姆霍兹方程在直角坐标系下的解

在直角坐标系下，亥姆霍兹方程表示为
$$
\left[ \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2} \right]p(x, y, z) + k^2 p(x, y, z) = 0
$$
采用分离变量法将声压的空间方程写为
$$
X''(x)Y(y)Z(z)+X(x)Y''(y)Z(z)+X(x)Y(y)Z''(z) + k^2X(x)Y(y)Z(z) = 0
$$
方程两边同时除以 $p(x,y,z) = X(x)Y(y)Z(z)$ 得
$$
\frac{X''(x)}{X(x)}+\frac{Y''(y)}{Y(y)} + \frac{Z''(z)}{Z(z)} + k^2 = 0
$$
令 $\frac{X''(x)}{X(x)}=-k_x, \frac{Y''(y)}{Y(y)}=-k_y, \frac{Z''(z)}{Z(z)}=-k_z$，且 $k_x^2 + k_y^2 + k_z^2 = -k^2$，可以推出
$$
X''(x)+k_x^2X(x) = 0 \\
Y''(y)+k_y^2Y(y) = 0 \\
Z''(z)+k_z^2Z(z) = 0 \\
$$
根据二阶常系数微分方程的解可得
$$
X(x) = Ae^{jk_xx} + Be^{-jk_xx} \\
Y(y) = Ce^{jk_yy} + De^{-jk_yy} \\
X(x) = Ee^{jk_zz} + Fe^{-jk_zz} \\
$$
则声场的空间解表示为
$$
p(x, y, z) = ACEe^{j(k_xx+k_yy+k_zz)}+ADEe^{j(k_xx-k_yy+k_zz)}+BCEe^{j(-k_xx+k_yy+k_zz)} \\
+ BDE e^{j(-k_xx-k_yy+k_zz)} + ACFe^{j(k_xx+k_yy-K_zz)} + ADFe^{j(k_xx-k_yy-K-zz)} \\
+BCFe^{j(-k_xx+k_yy-k_zz)}+BDFe^{j(-k_xx-k_yy-k_zz)} \\
 = C_1e^{j(k_xx+k_yy+k_zz)} + \cdots + C_8e^{j(-k_xx-k_yy-k_zz)} \\
 = \sum_{k_xx}\sum_{k_yy}C_{k_x, k_y}e^{\pm j \overrightarrow{k} \cdot \overrightarrow{r}}
$$
其中，$\overrightarrow{k} = (k_x, k_y, k_z), \overrightarrow{r} = (x, y, z)$。

再加上时间因子 $e^{j\omega t}$，直角坐标系中的波动方程解为
$$
p(\overrightarrow{r}, t) = \sum_{k_xx}\sum_{k_yy}C_{k_x, k_y}e^{j \left(\omega t\pm  \overrightarrow{k} \cdot \overrightarrow{r}\right)}
$$
直角坐标系下，波动方程形式解中的每一项，都表示向 $\overrightarrow{k}$ 方向传播的平面简谐波。

- **等相位面为平面；**
- **$\overrightarrow{k}$ 为等相位面法线方向；**
- **等相位面移动的速度为 $c_0$。**

<font color="red">任意空间分布的声波场（遵循亥姆霍兹方程），都可以分解为许多平面波的叠加形式——平面波分解。</font>

## 15.3 定解条件

**定解条件分为四类：边界条件、辐射条件、奇性条件和初始条件。**它们是声波传播过程中所满足的具体条件，通过这些条件，波动方程可求出特定解。

### 15.3.1 边界条件

1. <font color="red">**绝对软边界——声压释放边界**</font>

   - 第一类齐次性边界：$p = 0$；
   - 第一类非齐次性边界：$p = p_s$，当两种介质中的声压在分界面上连续时，$p_1 = p_2$。

2. <font color="red">**绝对硬边界——质点振速边界**</font>

   - 第二类齐次性边界：$\left.\overrightarrow{n} \cdot \overrightarrow{u} \right|_{z_0} = 0$；
   - 第二类非齐次性边界：$\left.\overrightarrow{n} \cdot \overrightarrow{u} \right|_{z_0} = u_s$，当两种介质在分界面处的法向振速连续时，$u_{1n} = u_{2n}$。

3. <font color="red">**混合边界条件——压力和质点振速线性组合**</font>
   $$
   \left. \frac{\partial p}{\partial n} + ap \right|_s = f(s)
   $$

   - 当 $a$ 为常数时，则称为第三类边界条件；
   - 若 $f(s) = 0$，则称为阻抗边界条件。

### 15.3.2 辐射条件

**当无穷远处没有声源存在时，声场应具有扩散波的性质，即无穷远处没有声场，声压为零。**

### 15.3.3 奇性条件

**均匀发散球面波在声源处存在奇异点，不满足波动方程。处理的方法是引入狄利克雷函数。**

### 15.3.4 初始条件

**当求远离初始时刻的稳态解时可不考虑初始条件。**
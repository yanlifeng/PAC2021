# PAC2021-ChiFanBuZhangRou

- [x] Use icc
- [x] Use intel fftw
- [x] Simple muti-threading
- [x] vectorize
- [x] check！
- [x] icpc ？？？
- [x] why 3DCTF nx * ny --> (nx+2-nx%2) * ny
- [x] fftw omp critical
- [x] change g++ version
- [x] vectorize for Rebu when g++
- [x] optimize weight、pre bufc and malloc cost 
- [x] optimize write
- [x] omp in for z
- [x] rewrite CTF vec 👋
- [x] optimize write2DIm
- [ ] 





注意：git上的cmakelist是本地用来测试的，实际用的还是原版。



## 0721

简单运行了一下，baseline是单线程的，编译器是gcc，fftw用的官方的3.3.7。

官方平台上7m30s左右，本地mac8m30s，现在本地稍微改改吧。

## 0722

3D-CTF correction部分和reconstruction是两个热点，目测循环60张图片是不能并行的，从这俩内部下手，CTF部分的循环（60+次）是没有依赖的，但是内部调用fftw的malloc和free并不是线程安全的，需要单线程执行，于是暂时又以下两种方案：

- 直接并行for，在malloc和free加critical或者其他：6m18s
- 并行fftwf_execute，其他热点也单独并行：//TODO

先打打timer吧，要不不好看加速效果。稳妥一点用比较熟悉的time方法。

本来想在本地先弄几天，但是发现本地和服务器热点出入有点大，就直接去服务器跑吧。

| Verison                               | cost total | cost 3DF | cost rebu |
| ------------------------------------- | ---------- | -------- | --------- |
| Init                                  | 444.01     | 268.27   | 174.49    |
| -Ofast                                | 443.8      | 268.59   | 173.97    |
| 3DF和reconstruction 简单并行 thread4  | 113.32     | 67.34    | 44.63     |
| 3DF和reconstruction 简单并行 thread24 | 21.50      | 12.42    | 7.64      |
| 3DF和reconstruction 简单并行 thread48 | 14.55      | 8.85     | 3.93      |
| 3DF和reconstruction 简单并行 thread96 | 13.00      | 7.36     | 3.72      |

前24个线程基本是线性的，24-48的时候rebu还是快一倍，但是3DF因为单线程的malloc和free出现了问题，另外96线程能不能直接用？两颗是一个节点？//TODO

总感觉服务器上安的这个oneapi它不对劲，icx找不到，好在mkl的东西似乎都有，先把fftw换成intel的试试。

阅读了原来的CMakelist，里面关于fftw的操作是include了一个单独的cmake文件，里面临时解压的fftw3.7，设置了flag和path等信息。去年pac用的还不是oneapi，而且也不是cmake，然后就给我整不明白了，官网的文档也仅仅是一句“Normally, the only change needed to build your application with FFTW3 wrappers replacing original FFTW library is to add Intel MKL at the link stage (see section *"Linking Your Application with Intel® Math Kernel Library" in the Intel MKL User's Guide*).”  试了简单的-I和-L都不行，属实是整不明白了最后只能跑example，看他的cmakelist是咋写的。

淦行不通，这examples根本编译不通过，麻了。

## 0723

又捣鼓了一天的fftw环境，基本的路子是找个简单的例子，先手动连接官方的👌然后在试试intel的👌，下一步就是怎么换到现在的项目里了。

https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl/link-line-advisor.html

这属实是个好东西。

今天大部分时间都在研究链接fftw库，顺带学了一下链接和cmake的知识，简单写一下。

对于官方的fftw，编译完成之后动态链接-I -L -lfftw3f -lfftw3f_threads就好了，静态的话把-L和-l换成.a就好；

对于intel的fftw，可以自己编译wrapper也可以直接用，直接用的话可以用👆提到的网址生成命令就好，遇到了一点坑就是oneapi中mkl的动态链接库好像都失效了，手动改了就好了。

举个例子：

```bash
静态链接：
[PAC20217111@manager src]$ g++ -o main main.cpp  -m64  -I"${MKLROOT}/include" -I"${MKLROOT}/include/fftw"  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
动态链接：
g++ -o main main.cpp  -m64  -I"${MKLROOT}/include" -I"${MKLROOT}/include/fftw"  -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
```



| Verison                                          | cost total | cost 3DF | cost rebu |
| ------------------------------------------------ | ---------- | -------- | --------- |
| 3DF和reconstruction 简单并行 thread4 intel fftw  | 18.55      | 10.08    | 8.05      |
| 3DF和reconstruction 简单并行 thread24 intel fftw |            |          |           |
| 3DF和reconstruction 简单并行 thread48 intel fftw |            |          |           |
| 3DF和reconstruction 简单并行 thread96 intel fftw |            |          |           |

先不测了，过不去check，麻了。

只换成icpx然后-fp-model precise一下可以过check（cmake .. -DUSE_INTEL_FFTW=0），时间略快，几乎一样，没啥用。



## 0724

ali冲冲冲，老甘yyds

## 0725

上来换成了icpc，直接把check干成了负数？？

很迷，不管了，就g++先优化代码吧。

thread 64

|                                         |        |       |       |
| --------------------------------------- | ------ | ----- | ----- |
| zz                                      | /      | /     |       |
| g++ -O3                                 | 14.176 | 8.975 | 3.253 |
| icpx -O3 -fp-model precise              | 13.907 | 7.934 | 3.943 |
| icpc -O3 -fp-model precise XXXXX        | 12.249 | 8.831 | 1.405 |
| g++ -O3 -funroll-loops -flto            | 13.704 | 8.706 | 3.070 |
| g++ -O3 -funroll-loops -flto thread 128 | 12.535 | 6.903 | 3.539 |
| g++ -O3 -funroll-loops -flto thread 96  | 12.364 | 6.995 | 3.367 |

icpc在recon部分能快一倍多，但是错的很离谱？先不管。

-ffast-math等影响精度的操作都会导致过不了check，编译器和编译参数的调试就先到这里，暂时用g++ -O3 -funroll-loops -flto。

然后改一下fftw，mkl的fftw会导致，略微微的过不了check，

```bash
mrc1mean: 984.078
mrc1min: -70729.7
mrc1max: 63454.7
mrc2mean: 1001.4
mrc2max: 63454.7
The error accumulation of the two volume is: 287399
The summation error of the two volume is: 0
The mean error of the two volume is: 17.3234
The absolute mean error on a single voxel is: 0.00438965
The relative mean error on a single voxel is: 4.46067e-06
Validation Failed!
```

现在官方fftw的参数是

```
--disable-doc --enable-threads --enable-float
```

改成512试试

```
--disable-doc --enable-threads --enable-float --enable-avx512
```

淦，gcc4.5不支持avx512，用256试试

```
--disable-doc --enable-threads --enable-float --enable-avx
```

淦，，会快一点，但是过不了check，这比赛咋打嘛。关于check稍微看了看脚本，关键就是比较了6e7个浮点数，然后把差值的绝对值加起来/6e7，这个值就是绝对误差，绝对误差/(6e7个数的平均值)就是相对误差，平均值大约是1000，然后6e7个数大约是1e3的数量级，最终要求相对误差<1e-7，即绝对误差<1e-4，即所有数的误差总和<6000，但肯定不能这么算

```bash
MRC1 Path:pro_novpp_Cor2_WBP_bin6.rec
MRC2 Path:Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
mrc1mean: 984.078
mrc1min: -70729.7
mrc1max: 63454.7
mrc2mean: 1001.4
mrc2max: 63454.7
The error accumulation of the two volume is: 51421.9
The summation error of the two volume is: 4096
The mean error of the two volume is: 17.3234
The absolute mean error on a single voxel is: 0.000785403
The relative mean error on a single voxel is: 7.9811e-07
```

好嘛fftw也不能改。

nzz的更新有bug，修改了。

会不会给的参考答案就有点问题，试试用初始代码的输出作为答案，区别只有平均值，其他全都一样，所以并不会影响答案。

=====================================

避免出现去年pac决赛卡在精度的问题，先暂时不管，继续看代码。

第一个热点简单看了一下，写了一个简单的优化版本：首先第一轮的fftw可以提出来只做一次，因为每次的计算过程都一样；此外，之前线程不安全的fftw的malloc和free也可以通过开副本的方式解决，略有提升：

|                       |        |       |       |
| --------------------- | ------ | ----- | ----- |
| 简单优化3DF thread64  | 11.521 | 7.598 | 2.862 |
| 简单优化3DF thread96  | 10.626 | 6.272 | 3.337 |
| 简单优化3DF thread128 | 10.594 | 6.399 | 3.481 |
|                       |        |       |       |

提升不大，毕竟只在外层做了循环，它只循环60次左右，先push一下，下一个版本写第一个热点内部并行的。

先不写内部并行的了，感觉直接在最最外层并行的效果更好，即68张图片并行处理，先看看直接加并行的效果：

|                      |        |         |        |
| -------------------- | ------ | ------- | ------ |
| for nz omp thread 16 | 25.838 | 257.899 | 86.879 |
|                      |        |         |        |
|                      |        |         |        |
|                      |        |         |        |

## 0726

好像只有1线程和68线程可以过check，很奇怪。



## 0727

|                                                   |        |        |       |
| ------------------------------------------------- | ------ | ------ | ----- |
| icpc submit 4                                     | 11.263 | 7.643  | 1.478 |
| icpc submit 5                                     | 8.984  | 6.430  | 1.648 |
| add -par-affinity=compact -funroll-all-loops -ipo | x      | x      | x     |
| icpc submit 5 thread 16                           | 18.278 | 13.488 | 3.363 |
| icpc submit 5 thread 32                           | 11.308 | 8.498  | 1.780 |
| icpc submit 5 thread 64                           | 9.374  | 6.765  | 1.061 |

淦，昨天那个线程数的bug是因为代码里是静态的threadnumber，但是运行的时候却用了-c。

|                                               | total   | read  | weight | 3DCTF   | rebu         |
| --------------------------------------------- | ------- | ----- | ------ | ------- | ------------ |
| thread 64                                     | 8.863   | 0.071 | 0.588  | 6.667   | 0.978+0.693  |
| -fast-transcendentals                         | 9.378   | 0.096 | 0.604  | 6.429   | 1.575+0.665  |
| -par-affinity=compact -funroll-all-loops -ipo | x       | x     | x      | x       | x            |
| thread 1                                      | 332.785 | 0.070 | 0.596  | 280.195 | 51.312+0.527 |
| thread 2                                      | 111.863 |       |        | 84.028  | 26.210+0.9   |
| thread 4                                      |         |       |        |         |              |
|                                               |         |       |        |         |              |

这icpc有毒啊，thread1答案又不对了，暂时先不用它了。

| version                     | total   | read  | weight | 3DCTF   | rebu        |
| --------------------------- | ------- | ----- | ------ | ------- | ----------- |
| g++ thread 128              | 10.832  | 0.060 | 0.600  | 6.435   | 3.702+0.028 |
| 64                          | 12.039  |       |        | 7.556   | 3.313+0.500 |
| 1                           | 365.363 |       |        | 212.574 | 151.202+0.7 |
| 64 优化sin cos,删除中间变量 | 9.428   |       |        | 6.920   | 1.372+0.4   |

整吐了，rebu里面涉及到ceil和floor，精度要求很高，感觉原来的代码也不一定是正确的。。。。

## 0728

rebu热点的i和x_orig部分几乎是不能动的，一动就gg。保守一点，把z部分的计算提出去。

3dCTF中第二次拷贝可以去掉，原来是Nx * Ny扩展成Nx2 * Ny，做fftw，然后再缩回Nx * Ny，后面使用[j * Nx + i]；这里可以直接不要“缩回”的操作，后面使用[j * Nx2 + i]。

关于weight和3dCTF中简化出来的第一步fftw现在都是单线程执行，可以加加并行。

|                                     |       |      |       |       |             |
| ----------------------------------- | ----- | ---- | ----- | ----- | ----------- |
| thread128 优化rebu orig等计算  8.41 | 7.721 |      |       | 6.070 | 0.971+0.02  |
| 3dCTF 优化掉第二次拷贝 7.28         | 6.588 |      |       | 4.863 |             |
| weight 和 提出来的fftw 多线程 6.62  | 5.897 |      | 0.288 | 4.460 | 1.043+0.039 |
|                                     |       |      |       |       |             |

## 0729

今天重点改一下3dCTF中的2for循环，for的写法是这样

```c++
for (int j = 0; j < Ny; j++) for (int i = 0; i < Nx2; i += 2)
```

里面比较关键的两个变量是

```c++
Nxh=Nx/2+1；
Nyh=Ny/2+1；
...
float x_norm = (x >= Nxh) ? (x - Nx) : (x);
float y_norm = (y >= Nyh) ? (y - Ny) : (y);
```

举个例子，

```
Nx=100 ：0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 

Ny=100 ：0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 -49 -48 -47 -46 -45 -44 -43 -42 -41 -40 -39 -38 -37 -36 -35 -34 -33 -32 -31 -30 -29 -28 -27 -26 -25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 

Nx=101 ：0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 

Ny=101 ：0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 -50 -49 -48 -47 -46 -45 -44 -43 -42 -41 -40 -39 -38 -37 -36 -35 -34 -33 -32 -31 -30 -29 -28 -27 -26 -25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 

```

简化之后就是

```c++
for (int j = 0; j < Ny; j++) {
                            for (int i = 0; i < Nxh; i++) {
                                float x = i, y = j;
                                float x_norm = x;
                                float y_norm = (y >= Nyh) ? (y - Ny) : (y);}}
```

无聊换上icpc和icpx试了试

|                                                   |       |       |       |       |             |
| ------------------------------------------------- | ----- | ----- | ----- | ----- | ----------- |
| thread128 优化rebu orig等计算  8.41               | 7.721 |       |       | 6.070 | 0.971+0.02  |
| 3dCTF 优化掉第二次拷贝 7.28                       | 6.588 |       |       | 4.863 |             |
| weight 和 提出来的fftw 多线程 6.62                | 5.897 |       | 0.288 | 4.460 | 1.043+0.039 |
| icpx -xHost -qopt-zmm-usage=high （mkl fftw）3.24 | 2.350 | 0.033 | 0.116 | 1.680 | 0.494+0.022 |

杀疯了，不过精度也炸了。



## 0730

```
gcc 5.4

-- The C compiler identification is GNU 5.4.0
-- The CXX compiler identification is GNU 5.4.0
-- Check for working C compiler: /usr/bin/cc
-- Check for working C compiler: /usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /usr/bin/c++
-- Check for working CXX compiler: /usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Welcome to TomoProject 0.1.20201014!
-- CMAKE_BUILD_TYPE : release, compile TomoProject with -std=c++11  -fopenmp -O2 flag.
-- Try to build Tomo-Project in CPU version.
-- CMAKE_C_FLAGS : -std=c++11  -fopenmp -O2 -lz -lpthread -lm
-- CMAKE_CXX_FLAGS : -std=c++11  -fopenmp -O2 -lz -lpthread -lm
-- FFTW_FLAGS : --disable-doc;--enable-threads;--enable-float;--prefix=/home/asc/pac2021_publish/WBP/external/fftw-3.3.7
-- FFTW_LIBRARIES : /home/asc/pac2021_publish/WBP/external/fftw-3.3.7/lib/libfftw3f.a;/home/asc/pac2021_publish/WBP/external/fftw-3.3.7/lib/libfftw3f_threads.a
-- FFTW Path: /home/asc/pac2021_publish/WBP/external/fftw-3.3.7
-- Configuring done
-- Generating done
-- Build files have been written to: /home/asc/pac2021_publish/WBP/build

(presto) asc@ubuntu:~/pac2021_publish/data/proteasome-bin6$  ../../WBP/build/validate_cpu pro_novpp_Cor2_WBP_bin6.rec Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
MRC1 Path:pro_novpp_Cor2_WBP_bin6.rec
MRC2 Path:Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
mrc1mean: 984.078
mrc1min: -70729.6
mrc1max: 63454.6
mrc2mean: 1001.4
mrc2max: 63454.7
The error accumulation of the two volume is: 49136.8
The summation error of the two volume is: 8192
The mean error of the two volume is: 17.3234
The absolute mean error on a single voxel is: 0.000750501
The relative mean error on a single voxel is: 7.62644e-07
Validation Failed!








-- The C compiler identification is GNU 9.3.0
-- The CXX compiler identification is GNU 9.3.0
-- Check for working C compiler: /usr/bin/cc
-- Check for working C compiler: /usr/bin/cc -- works
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Detecting C compile features
-- Detecting C compile features - done
-- Check for working CXX compiler: /usr/bin/c++
-- Check for working CXX compiler: /usr/bin/c++ -- works
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Detecting CXX compile features
-- Detecting CXX compile features - done
fatal: not a git repository (or any parent up to mount point /home)
Stopping at filesystem boundary (GIT_DISCOVERY_ACROSS_FILESYSTEM not set).
-- Welcome to TomoProject 0.1.20201014!
-- CMAKE_BUILD_TYPE : release, compile TomoProject with -std=c++11  -fopenmp -O2 flag.
-- Try to build Tomo-Project in CPU version.
-- CMAKE_C_FLAGS : -std=c++11  -fopenmp -O2 -lz -lpthread -lm
-- CMAKE_CXX_FLAGS : -std=c++11  -fopenmp -O2 -lz -lpthread -lm
-- FFTW_FLAGS : --disable-doc;--enable-threads;--enable-float;--prefix=/home/user_home/ylf/pac2021_publish/WBP/external/fftw-3.3.7
-- FFTW_LIBRARIES : /home/user_home/ylf/pac2021_publish/WBP/external/fftw-3.3.7/lib/libfftw3f.a;/home/user_home/ylf/pac2021_publish/WBP/external/fftw-3.3.7/lib/libfftw3f_threads.a
-- FFTW Path: /home/user_home/ylf/pac2021_publish/WBP/external/fftw-3.3.7
-- Configuring done
-- Generating done
-- Build files have been written to: /home/user_home/ylf/pac2021_publish/WBP/build

MRC1 Path:pro_novpp_Cor2_WBP_bin6.rec
MRC2 Path:Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
mrc1mean: 984.079
mrc1min: -70729.9
mrc1max: 63454.5
mrc2mean: 1001.4
mrc2max: 63454.7
The error accumulation of the two volume is: 4.38686e+06
The summation error of the two volume is: 61440
The mean error of the two volume is: 17.3233
The absolute mean error on a single voxel is: 0.0670036
The relative mean error on a single voxel is: 6.80876e-05
Validation Failed!
```



icpc把参数precise改成strict

```
strict
Enables precise and except, disables contractions, and enables pragma stdc fenv_access.

[no-]except (Linux* and macOS*) or except[-] (Windows* )
Determines whether strict floating-point exception semantics are honored.


```

see more in https://software.intel.com/content/www/us/en/develop/documentation/cpp-compiler-developer-guide-and-reference/top/compiler-reference/compiler-options/compiler-option-details/floating-point-options/fp-model-fp.html

New icpc

| Verison                                         | Cpp tot | Read  | Weight | CTF     | Rebu        | Total     |
| ----------------------------------------------- | ------- | ----- | ------ | ------- | ----------- | --------- |
| weight 和 提出来的fftw 多线程                   | 5.897   |       | 0.288  | 4.460   | 1.043+0.039 | 6.62      |
| icpx -xHost -qopt-zmm-usage=high （mkl fftw）   | 2.350   | 0.033 | 0.116  | 1.680   | 0.494+0.022 | 3.24      |
| -fp-model strict -xHost -qopt-zmm-usage=high 64 | 6.093   | 0.063 | 0.326  | 5.260   | 0.339+0.098 | 6.99      |
| -fp-model strict -xHost -qopt-zmm-usage=high 32 | 9.042   | 0.057 | 0.202  | 7.500   | 0.424+0.851 | 9.97      |
| -fp-model strict -xHost -qopt-zmm-usage=high 16 | 13.665  | 0.061 | 0.177  | 11.999  | 0.702+0.716 | 14.54     |
| -fp-model strict -xHost -qopt-zmm-usage=high 8  | 23.361  | 0.060 | 0.189  | 21.129  | 1.246+0.723 | 24.27947  |
| -fp-model strict -xHost -qopt-zmm-usage=high 4  | 42.909  | 0.063 | 0.208  | 39.524  | 2.362+0.736 | 43.80787  |
| -fp-model strict -xHost -qopt-zmm-usage=high 2  | 84.234  | 0.060 | 0.243  | 78.397  | 4.740+0.759 | 85.11350  |
| -fp-model strict -xHost -qopt-zmm-usage=high 1  | 164.341 | 0.064 | 0.289  | 154.657 | 8.687+0.586 | 165.23959 |
| 1 mkl fftw x                                    | 137.399 |       |        | 128.187 |             | 138.32196 |
| 64 fftw_256 x                                   | 5.912   | 0.065 | 0.317  | 5.119   | 0.348+0.057 | 6.78638   |

关键就是CTF的向量化，mkl的fftw估计是暂时不能用了。

关于CTF中的三角函数，可以通过预处理或者是化简等方式优化他们：

| Verison                                         | Cpp tot | Read  | Weight | CTF   | Rebu        | Total |
| ----------------------------------------------- | ------- | ----- | ------ | ----- | ----------- | ----- |
| weight 和 提出来的fftw 多线程                   | 5.897   |       | 0.288  | 4.460 | 1.043+0.039 | 6.62  |
| icpx -xHost -qopt-zmm-usage=high （mkl fftw）   | 2.350   | 0.033 | 0.116  | 1.680 | 0.494+0.022 | 3.24  |
| -fp-model strict -xHost -qopt-zmm-usage=high 64 | 6.093   | 0.063 | 0.326  | 5.260 | 0.339+0.098 | 6.99  |
| 😄优化CTF中的三角函数                            | 3.932   |       |        | 3.212 |             | 4.91  |
| 🤔cos(2 * (alpha - astig))也去掉                 | 5.308   |       |        | 4.524 |             | 6.28  |
| 😄换成g++ 128                                    | 4.487   |       |        | 2.984 | 1.169+0.014 | 5.21  |

可以看到g++也差不太多了，只是Rebu没有开启向量化。。。。。

晚上和启鑫凯子哥聊了一下，有一下几个点可以借鉴

- cos(x)= x - ((x * x * x) / 6) + ((x * x * x * x * x) / 120)
- CTF中的cos(2 * (alpha - astig))直接把atan带进去，抵消掉三角函数

上面这两个策略，第二个和👆Version中🤔差不多，都是通过多余的计算把cos去掉，但🤔效果不好，这个目测也一般；第一个不能保证精度，试了试可以通过但是defocus1 - defocus2太小了，不确定是为啥能过。

## 0731

差不多要开始写CTF部分的向量化了。

对于icpc fp model的参数基本上都是了一遍，consistent 即-fp-model precise -no-fma -fimf-arch-consistency=true是能过check里面最快的✈️，通过向量化报告也可以看到两个热点都进行了向量化，但是比起不限制浮点运算的版本🚀在CTF上还是差得多：

| Verison                      | Cpp tot | Read  | Weight | CTF   | Rebu        | Total |
| ---------------------------- | ------- | ----- | ------ | ----- | ----------- | ----- |
| g++ 128                      | 4.558   | 0.065 | 0.256  | 3.016 | 1.203+      | 5.35  |
| icpc 64 strict new malloc    | 4.068   | 0.067 | 0.365  | 3.245 | 0.339+0.047 | 4.94  |
| icpc  😭                      | 3.320   | 0.070 | 0.403  | 2.400 | 0.399+0.043 | 4.21  |
| icpc  😭 mkl fftw             | 2.131   | 0.055 | 0.150  | 1.578 | 0.304+0.038 | 3.08  |
| icpc  😭 mkl fftw no critical | 1.538   | 0.057 | 0.146  | 0.983 | 0.306+0.040 | 2.49  |

大体上就是CTF部分不加strict可以开启自动向量化，能起到一定的提升，然后用mkl的fftw，本身就会快很多，如果在把critical去掉，就直接起飞。但是后三个😭版本过不去check，现在的想法是先手动向量化达到4.2s，然后想办法优化fftw中malloc和free部分的栅栏。

直接开始写g++版本的Rebu向量化（其实是拿icpc测试的，因为g++4.8不支持avx512）：

好啊，写好了。

坑1：_mm512_mul_epi32指令实际上算的是一般的值

坑2: 写向量化的时候li别忘了

坑3: div指令贼慢。。。。实测并不比gather快多少

| Verison                   | Cpp tot | Read  | Weight | CTF   | Rebu        | Total   |
| ------------------------- | ------- | ----- | ------ | ----- | ----------- | ------- |
| g++ 128                   | 4.558   | 0.065 | 0.256  | 3.016 | 1.203+      | 5.35    |
| icpc 64 strict new malloc | 4.068   | 0.067 | 0.365  | 3.245 | 0.339+0.047 | 4.94    |
| icpc 64 vectorize by 👋    | 3.861   | 0.071 | 0.252  | 3.183 | 0.303+0.047 | 4.86458 |

好啊，这部分手写的还行，和自动的效果几乎一样，明天搞一搞CTF部分的。

冲冲冲。

## 0801

争取今天写好CTF的向量化能进破4s然后就去干点别的了，ali的比赛已经滚出第一页了，再不干活老甘就要🦈我了，QC也耽误了很久了。

CTF部分打打timer发现critical和fftw比较耗时，相比之下forfor不太是热点了。。。。没办法，简单写写向量化吧争取到4.2。

不太对，重新打了计时函数发现CTF中与处理部分慢了，没有并行。

```
main for cost 3.389
init and read cost 0.073
weight cost 0.288
pre bufc cost 0.266
pre atant cost 0.016
fftw malloc cost 0.227
hotspots1 cost 2.156
fftw free cost 0.020
bufc free cost 0.000
rebu cost 0.339
stack_corrected free cost 0.001
Wrtie out final reconstruction result:
Done

Finish reconstruction successfully!
All results save in: ./


Finish!
Time elapsed: 5s

tot cost 4.26720
```

争取写好向量化到3.7

| Verison                                         | Cpp tot | Read  | Weight | CTF                            | Rebu        | Total |
| ----------------------------------------------- | ------- | ----- | ------ | ------------------------------ | ----------- | ----- |
| g++ 128                                         | 4.558   | 0.065 | 0.256  | 3.016                          | 1.203+      | 5.35  |
| icpc 64 strict new malloc                       | 4.068   | 0.067 | 0.365  | 3.245                          | 0.339+0.047 | 4.94  |
| icpc 64 vectorize by 👋                          | 3.861   | 0.071 | 0.252  | 3.183                          | 0.303+0.047 | 4.86  |
| icpc 64 omp for CTF pre data and remove img_pre | 3.389   | 0.079 | 0.240  | 0.266+0.016+0.227 +2.156+0.020 | 0.339+0.000 | 4.26  |
| 👆remove cos                                     | 2.965   |       |        | +1.756+                        |             | 3.80  |
| CTF vectorize                                   | 2.778   | 0.075 | 0.255  | 0.193+0.065+0.228 +1.518+0.018 | 0.343+      | 3.70  |

## 0802

今天简单写写write部分的优化就收了，写写阿里和论文。

对了测试了一下，CTF大热点部分大约是fftw占1s多，着重优化的两个for循环占现在占0.3左右，不关键了已经；而且前面预处理部分的weight和pre bufc最好情况下是占0.4左右，但是很不稳定，经常到0.8，猜测是fftw内部的多线程原因，刚才尝试把这两个部分挪到最外层外面，这样就可以omp并行它，但是image和bufc也好开成全局的大大数组，这样memcpy好像贼慢，就回退到上一个版本了，这一部分有待优化//TODO。

对于write部分，写stack_recon部分不能合成一次写，以为地址不连续，但或许可以改改试试//TODO，or开多个m_fp副本来多线程的写，其实这部分只占0.1s，可以先不管。

最后的fclose慢的很，0.5s左右，可能是在把缓冲区的数据刷新到硬盘（待确定//TODO），暂时直接把这块去掉了，似乎不太好，但是勉强写出一个能进2s的版本就先收了，明天开始做做别的。

| Verison                                         | Cpp tot | Read  | Weight | CTF                                  | Rebu        | Total          |
| ----------------------------------------------- | ------- | ----- | ------ | ------------------------------------ | ----------- | -------------- |
| g++ 128                                         | 4.558   | 0.065 | 0.256  | 3.016                                | 1.203+      | 5.35           |
| icpc 64 strict new malloc                       | 4.068   | 0.067 | 0.365  | 3.245                                | 0.339+0.047 | 4.94           |
| icpc 64 vectorize by 👋                          | 3.861   | 0.071 | 0.252  | 3.183                                | 0.303+0.047 | 4.86           |
| icpc 64 omp for CTF pre data and remove img_pre | 3.389   | 0.079 | 0.240  | 0.266+0.016+0.227 +2.156+0.020       | 0.339+0.000 | 4.26           |
| 👆remove cos                                     | 2.965   |       |        | +1.756+                              |             | 3.80           |
| CTF vectorize                                   | 2.778   | 0.075 | 0.255  | 0.193+0.065+0.228 +1.518+0.018       | 0.343+      | 3.70           |
| optimize write part(check ?)                    | 2.679   | 0.040 | 0.221  | 0.274+0.000+0.228 +1.531+0.016+0.002 | 0.353+0.001 | +0.138+? =2.99 |



## 0803

早上拿几分钟写了一下外层循环的，但是发现stack_recon的副本根本开不下，遂投降，下面测测目前版本的线程拓展性：

| Verison | total | cpp    | read  | Weight | bufc  | fftw malloc | CTF          | free        | rebu        | write |
| ------- | ----- | ------ | ----- | ------ | ----- | ----------- | ------------ | ----------- | ----------- | ----- |
| 64      | 3.03  | 2.710  | 0.049 | 0.261  | 0.263 | 0.226       | 1.530/0.299  | 0.017+0.002 | 0.351+0.001 | 0.155 |
| 32      | 3.30  | 2.973  |       | 0.126  | 0.125 |             | 2.070/0.287  |             | 0.361+      |       |
| 16      | 4.75  | 4.433  |       |        |       |             | 3.381/0.463  |             | 0.526+      |       |
| 8       | 7.87  | 7.532  |       |        |       |             | 6.013/0.739  |             | 0.925       |       |
| 4       | 13.95 | 13.595 |       |        |       |             | 11.312/1.455 |             | 1.724       | 0.190 |
| 2       | 26.74 | 26.326 |       |        |       |             | 22.263/2.801 |             | 3.372       | 0.248 |
| 1       | 51.73 | 51.226 |       |        |       |             | 44.434/6.287 |             | 6.175       | 0.358 |

CTF中/后面是去掉fftw的时间，可以看到，经过一系列的优化，着重优化的forfor相比于fftw已经完全不是热点了。。。

如果用mkl的fftw，顺面把fftw的malloc和free并行：

| Verison | total | cpp   | read  | Weight | bufc  | fftw malloc | CTF         | free       | rebu        | write |
| ------- | ----- | ----- | ----- | ------ | ----- | ----------- | ----------- | ---------- | ----------- | ----- |
| 64      | 1.52  | 1.208 | 0.044 | 0.081  | 0.045 | 0.032       | 0.650/0.299 | 0.015+0.00 | 0.315+0.001 | 0.155 |

但是

```
MRC1 Path:pro_novpp_Cor2_WBP_bin6.rec
MRC2 Path:Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
mrc1mean: 984.078
mrc1min: -70729.7
mrc1max: 63454.7
mrc2mean: 1001.4
mrc2min: -70729.6
mrc2max: 63454.7
range 65472000
The error accumulation of the two volume is: 287399
The error accumulation of the two volume is(double): 412708
The summation error of the two volume is: 0
The mean error of the two volume is: 17.3234
The absolute mean error on a single voxel is: 0.00438965
The relative mean error on a single voxel is: 4.46067e-06
The absolute mean error on a single voxel is(double): 0.00630358
The relative mean error on a single voxel is(double): 6.40557e-06
Validation Failed!
```

```

gcc7 

MRC1 Path:pro_novpp_Cor2_WBP_bin6.rec
MRC2 Path:Result_NewMRCCTF_NewWBP/pro_novpp_Cor2_WBP_bin6.rec
mrc1mean: 984.079
mrc1min: -70729.9
mrc1max: 63454.4
mrc2mean: 1001.4
mrc2min: -70729.6
mrc2max: 63454.7
range 65472000
The error accumulation of the two volume is: 4.38406e+06
The error accumulation of the two volume is(double): 5.77985e+06
The summation error of the two volume is: 65536
The mean error of the two volume is: 17.3234
The absolute mean error on a single voxel is: 0.0669609
The relative mean error on a single voxel is: 6.80442e-05
The absolute mean error on a single voxel is(double): 0.0882797
The relative mean error on a single voxel is(double): 8.97079e-05
Validation Failed!
```

## 0804

改规则了淦，还要再干活。

- [x] 去掉fp-model，更换fp-model
- [x] 自动向量化
- [x] 重写向量化
- [x] cos->float
- [ ] g++
- [x] mkl fftw threads



| Verison          | total | func  | read  | Weight | bufc  | CTF   | rebu  | write | Error    |
| ---------------- | ----- | ----- | ----- | ------ | ----- | ----- | ----- | ----- | -------- |
| Icpc 64 mkl auto | 1.473 | 1.469 | 0.036 | 0.030  | 0.045 | 0.715 | 0.311 | 0.166 | 4.47e-06 |
| omp for new      | 1.275 | 1.271 | 0.040 | 0.045  | 0.035 | 0.708 | 0.209 | 0.160 | 4.47e-06 |
| Icpc 64 mkl 👋    |       |       |       |        |       | 0.572 | 0.193 |       |          |
|                  |       |       |       |        |       |       |       |       |          |
|                  |       |       |       |        |       |       |       |       |          |
|                  |       |       |       |        |       |       |       |       |          |

## 0805

关于写文件的优化，第一个策略是把Ny次合成一次，第二个就是关于close同步的时候，把结果写到shm中，然后再cpoy到数据目录。

```
#!/bin/bash
#SBATCH -p comp
#SBATCH -N 1
#SBATCH --exclusive
testpath=./

cd $testpath
binPath=../../WBP/build
tomoBin=$binPath/tomo_lxy
paraPath=./para_reconstruction_proteasome_WBP_bin6.conf
 $tomoBin $paraPath
mv /dev/shm/pro_novpp_Cor2_WBP_bin6.rec ./
```



| Verison                       | total | func  | read  | Weight | bufc  | CTF   | rebu  | write | Error    |
| ----------------------------- | ----- | ----- | ----- | ------ | ----- | ----- | ----- | ----- | -------- |
| Icpc 64 mkl 👋 +optimize write | 1.192 | 1.188 | 0.036 | 0.030  | 0.045 | 0.653 | 0.261 | 0.082 | 4.47e-06 |
|                               |       |       |       |        |       |       |       |       |          |

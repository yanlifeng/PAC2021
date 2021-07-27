# PAC2021-ChiFanBuZhangRou

- [x] Use icc
- [x] Use intel fftw
- [x] Simple muti-threading
- [ ] vectorize
- [ ] check！
- [ ] icpc ？？？



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
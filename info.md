# PAC2021-ChiFanBuZhangRou

- [x] Use icc
- [x] Use intel fftw
- [x] Simple muti-threading
- [ ] vectorize
- [ ] check！



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
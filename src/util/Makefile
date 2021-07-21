#笔记：include关键字：引入别的Makefile 可以引入多个Makefile，可以使用通配符
include ../../Makefile.inc

#patsubst笔记：
#	$(patsubst <pattern>,<replacement>,<text>)
#	名称：模式字符串替换函数。
#	功能：查找 <text> 中的单词（单词以“空格”、“Tab”或“回车”“换行”分隔）是否符合模式 <pattern> ，
#	如果匹配的话，则以 <replacement> 替换。这里， <pattern> 可以包括通配符 % ，表示任意长度的字串。如果 <replacement> 中也包含 % ，那么， <replacement> 中的这个 % 将是 <pattern> 中的那个 % 所代表的字串。（可以用 \ 来转义，以 \% 来表示真实含义的 % 字符）
#	返回：函数返回被替换过后的字符串。

#wildcard笔记：
#在Makefile规则中，通配符会被自动展开。但在变量的定义和函数引用时，通配符将失效。这种情况下如果需要通配符有效，就需要使用函数“wildcard”，它的用法是：$(wildcard PATTERN...) 。

UTIL_MODULE := $(patsubst %.c, %_c.o, $(wildcard *.c))
UTIL_MODULE += $(patsubst %.cc, %_cc.o, $(wildcard *.cc))
UTIL_MODULE += $(patsubst %.cpp, %_cpp.o, $(wildcard *.cpp))
UTIL_MODULE += $(patsubst %.cxx, %_cxx.o, $(wildcard *.cxx))

#subst笔记
#	$(subst <from>,<to>,<text>)
#	功能：把字串 <text> 中的 <from> 字符串替换成 <to> 。
UTIL_PROJECT_ROOT_DIR := $(subst src/util,,$(CURRENT_DIR))

#这个很重要，否则每次会因为找不到相应的文件而重新编译
#笔记：Makefile文件中的特殊变量 VPATH 就是完成这个功能的，如果没有指明这个变量，make只会在当前的目录中去找寻依赖文件和目标文件。
#	   如果定义了这个变量，那么，make就会在当前目录找不到的情况下，到所指定的目录中去找寻文件了。
VPATH = ../../obj/:../../lib
#下面这两行与上面这行等价
#vpath %.o $(UTIL_PROJECT_ROOT_DIR)/obj/
#vpath %.a $(UTIL_PROJECT_ROOT_DIR)/lib/
#笔记： 
#	vpath <pattern> <directories> 
#	vapth使用方法中的<pattern>需要包含 % 字符。 % 的意思是匹配零或若干字符，（需引用 % ，使用 \ ）
#	例如， %.h 表示所有以 .h 结尾的文件。
#	<pattern>指定了要搜索的文件集，
#	<directories>则指定了< pattern>的文件集的搜索的目录。例如：

#赋值运算符笔记：
#= 是最基本的赋值
#:= 是覆盖之前的值
#?= 是如果没有被赋值过就赋予等号后面的值
#+= 是添加等号后面的值
#Reference:https://www.cnblogs.com/wanqieddy/archive/2011/09/21/2184257.html

#addprefix笔记：
#$(addprefix <prefix>,<names...>)
#名称：加后缀函数——addsuffix。
#功能：把后缀 <suffix> 加到 <names> 中的每个单词后面。
UTIL_MODULE_TO_ARCHIVE := $(addprefix $(UTIL_PROJECT_ROOT_DIR)/obj/, $(UTIL_MODULE))


STATIC_LIB_NAME := libutil.a
#重要 编译命令 包括include path，第三方库的include path，library path
define compile-cmd
	$(CC) $(CFLAGS) -I$(INC_PATH) -I$(EXTERNAL_INC_PATH) -I$(EXTERNAL_INC_EIGEN_PATH) -L$(LIB_PATH) -L$(EXTERNAL_LIB_PATH) -o $(addprefix ../../obj/, $@) -c $<
endef
#笔记 $@, $<, $*等等都是自动化变量

#重要 库生成命令
define archive-cmd
	#$(AR) $(ARFLAGS) -o $(addprefix ../../lib/, $@) $(addprefix ../../obj/, $^)
	$(AR) $(ARFLAGS) -o $(addprefix ../../lib/, $@) $(UTIL_MODULE_TO_ARCHIVE) 
endef
all:$(UTIL_MODULE) $(STATIC_LIB_NAME)
%_c.o:%.c
	$(compile-cmd)
%_cc.o:%.cc
	$(compile-cmd)
%_cpp.o:%.cpp
	$(compile-cmd)
%_cxx.o:%.cxx
	$(compile-cmd)
$(STATIC_LIB_NAME):$(UTIL_MODULE)
	$(archive-cmd)
.phony:clean
clean:
	rm -rf ../../lib/$(STATIC_LIB_NAME)


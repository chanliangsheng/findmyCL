QT       += core gui
QT += sql

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    cardiolipin/CL/src/CL.cpp \
    cardiolipin/CL/src/ClSpecificStructure.cpp \
    cardiolipin/DLCL/src/DLCL.cpp \
    cardiolipin/DLCL/src/DlclSpecificStructure.cpp \
    cardiolipin/MLCL/src/MLCL.cpp \
    cardiolipin/MLCL/src/MlclSpecificStructure.cpp \
    cardiolipin/src/PaNode.cpp \
    cardiolipin/src/cardiolipin.cpp \
    database/src/database.cpp \
    database/src/databaserecord.cpp \
    main.cpp \
    mainwindow.cpp \
    mzml/base64/src/base64.cpp \
    mzml/ms1/src/ms1.cpp \
    mzml/ms2/src/headgroup.cpp \
    mzml/ms2/src/ms2.cpp \
    mzml/ms2/src/pa.cpp \
    mzml/src/mzml.cpp \
    workflow/FragmentCombiner/src/FragmentCombiner.cpp \
    workflow/FragmentFinder/src/FragmentFinder.cpp \
    workflow/HeadgroupFinder/src/HeadgroupFinder.cpp \
    workflow/Ms1LibraryMatcher/src/Ms1LibraryMatcher.cpp \
    workflow/MsLevelMatcher/src/MsLevelMatcher.cpp \
    workflow/src/workflow.cpp

HEADERS += \
    cardiolipin/CL/include/CL.h \
    cardiolipin/CL/include/ClSpecificStructure.h \
    cardiolipin/DLCL/include/DLCL.h \
    cardiolipin/DLCL/include/DlclSpecificStructure.h \
    cardiolipin/MLCL/include/MLCL.h \
    cardiolipin/MLCL/include/MlclSpecificStructure.h \
    cardiolipin/include/PaNode.h \
    cardiolipin/include/cardiolipin.h \
    database/include/database.h \
    database/include/databaserecord.h \
    mainwindow.h \
    mzml/base64/include/base64.h \
    mzml/include/mzml.h \
    mzml/ms1/include/ms1.h \
    mzml/ms2/include/headgroup.h \
    mzml/ms2/include/ms2.h \
    mzml/ms2/include/pa.h \
    workflow/FragmentCombiner/include/FragmentCombiner.h \
    workflow/FragmentFinder/include/FragmentFinder.h \
    workflow/HeadgroupFinder/include/HeadgroupFinder.h \
    workflow/Ms1LibraryMatcher/include/Ms1LibraryMatcher.h \
    workflow/MsLevelMatcher/include/MsLevelMatcher.h \
    workflow/include/workflow.h

FORMS += \
    mainwindow.ui

INCLUDEPATH += $$PWD/mzml/include
INCLUDEPATH += $$PWD/mzml/ms1/include
INCLUDEPATH += $$PWD/mzml/ms2/include
INCLUDEPATH += $$PWD/mzml/base64/include
INCLUDEPATH += $$PWD/database/include
INCLUDEPATH += $$PWD/cardiolipin/include
INCLUDEPATH += $$PWD/cardiolipin/CL/include
INCLUDEPATH += $$PWD/cardiolipin/MLCL/include
INCLUDEPATH += $$PWD/cardiolipin/DLCL/include
INCLUDEPATH += $$PWD/workflow/include
INCLUDEPATH += $$PWD/workflow/Ms1LibraryMatcher/include
INCLUDEPATH += $$PWD/workflow/MsLevelMatcher/include
INCLUDEPATH += $$PWD/workflow/HeadgroupFinder/include
INCLUDEPATH += $$PWD/workflow/FragmentFinder/include
INCLUDEPATH += $$PWD/workflow/FragmentCombiner/include

# 导入第三方库的头文件和动态库地址
INCLUDEPATH += 3rdparty/MyLibary/include
LIBS += -L $$PWD/3rdparty/MyLibary/lib -l tinyxml2
LIBS += -lz

# 把库文件加入到运行时的文件夹中
MY_PWD = $$PWD
OUTPUT_PWD = $$OUT_PWD
# 将变量MY_PWD中的正斜杠替换为反斜杠
MY_PWD_WIN = $$replace(MY_PWD, /, \\)
OUTPUT_PWD_WIN = $$replace(OUTPUT_PWD, /, \\)

QMAKE_PRE_LINK = xcopy $$MY_PWD_WIN\3rdparty\MyLibary\lib\tinyxml2.dll $$OUTPUT_PWD_WIN\debug /y

# 复制db文件夹到编译后项目中
MY_DB_PWD = $$MY_PWD_WIN\database.db
OUT_DB_PWD = $$OUTPUT_PWD_WIN\

copy_db = $$system("echo d| xcopy $${MY_DB_PWD} $${OUT_DB_PWD}")


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

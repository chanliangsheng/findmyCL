#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDebug>
#include <base64.h>
#include <database.h>
#include <mzml.h>
#include <Ms1LibraryMatcher.h>
#include <database.h>
#include <MsLevelMatcher.h>
#include <HeadgroupFinder.h>
#include <FragmentFinder.h>
#include <FragmentCombiner.h>
#include <QFileDialog>
#include <QThread>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

public:
    QThread *thread_1;

//信号
signals:
    void SendMs1CsvFileName(QString file_name);
    void SendMs2MzmlFileNames(QStringList file_names);
};
#endif // MAINWINDOW_H

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

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

#include "mainwindow.h"
#include "ui_mainwindow.h"



using namespace std;

//// dda
//MainWindow::MainWindow(QWidget *parent)
//    : QMainWindow(parent)
//    , ui(new Ui::MainWindow)
//{
//    ui->setupUi(this);

//    setWindowTitle(tr("心磷脂定性软件")); // 设置窗口标题为"My Widget"

//    thread_1 = new QThread(this);//新建一个线程
//    //新建mzml对象，移动到线程thread_1，线程启动
//    Mzml* mzml_ptr = new Mzml();
//    Database* database_ptr = new Database();

//    Ms1LibraryMatcher* ms1_library_matcher_ptr = new Ms1LibraryMatcher();

//    MsLevelMatcher* ms_level_matcher_ptr = new MsLevelMatcher();
//    HeadgroupFinder* headgroup_finder_ptr = new HeadgroupFinder();
//    FragmentFinder* fragment_finder_ptr = new FragmentFinder();
//    FragmentCombiner* fragment_combiner_ptr = new FragmentCombiner();
//    mzml_ptr->moveToThread(thread_1);
//    ms1_library_matcher_ptr->moveToThread(thread_1);
//    database_ptr->moveToThread(thread_1);
//    ms_level_matcher_ptr->moveToThread(thread_1);
//    headgroup_finder_ptr->moveToThread(thread_1);
//    fragment_finder_ptr->moveToThread(thread_1);
//    fragment_combiner_ptr->moveToThread(thread_1);
//   //ui

//    this->ui->TableWidgetMatchMs1Result->horizontalHeader()->setStretchLastSection(true);
//    this->ui->TableWidgetMatchMs1Result->horizontalHeader()->setVisible(true);
//    this->ui->TableWidgetSpliceResult->horizontalHeader()->setVisible(true);

//    //打开csv文件，发射信号
//    connect(this->ui->ActionOpenMs1 , &QAction::triggered , this , [=](){
//        QFileDialog dialog(this);
//        dialog.setDirectory("D:/test");
//        QString file_name = dialog.getOpenFileName();//获取文件路径
//        emit this->SendMs1CsvFileName(file_name);
//    });
//    //打开多个mzml文件，发射信号
//    connect(this->ui->ActionOpenMs2 , &QAction::triggered , this , [=]{
//        QFileDialog dialog(this);
//        dialog.setDirectory("D:/findmyCL/data/profile");
//        QStringList file_names =  dialog.getOpenFileNames();
//        emit this->SendMs2MzmlFileNames(file_names);
//    });



//    //读取一级和二级
//    //读取一级csv文件
//    connect(this,&MainWindow::SendMs1CsvFileName , mzml_ptr,[=](QString file_name){
//        mzml_ptr->ReadMs1FromCsv(file_name);

////        //idx
////        mzml_ptr->ConvertMs1RtUnit("mintosec");

//        QString message = "Ms1 count:" + QString::number(mzml_ptr->GetLocalMs1Vector().size());
//        this->ui->LabelMs1Message->setText(message);

//    });
//    //读取二级mzml文件，加载数据库
//    connect(this,&MainWindow::SendMs2MzmlFileNames , mzml_ptr,[=](QStringList file_names){
//        mzml_ptr->ReadMs2FromMzmls(file_names);

////        mzml_ptr->ConvertMs2RtUnit("mintosec");//转换为秒

//        database_ptr->LoadAllTable();
//        QString message = "Ms2 count:" + QString::number(mzml_ptr->GetLocalMs2Vector().size());
//        this->ui->LabelMs2Message->setText(message);

//        mzml_ptr->DeleteMs2LowIntensityFragment(0.1);
//    });


//    //一级配对
//    connect(this->ui->PushButtonMatchMs1 , &QPushButton::clicked , mzml_ptr , [=]() mutable {
//        ms1_library_matcher_ptr->m_ppm = this->ui->DoubleSpinBoxMatchMs1PPm->value();
//        ms1_library_matcher_ptr->m_tolerance_rt_scope = this->ui->DoubleSpinBoxMatchMs1ToleranceRt->value();
//        ms1_library_matcher_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxMatchMs1PPmWithHalfScore->value();
//        ms1_library_matcher_ptr->MatchMs1WithAllTables(*mzml_ptr , *database_ptr);
//        ms1_library_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1Result);//打印结果

//    });
//    //根据Rt打印结果
//    connect(this->ui->PushButtonMatchMs1FilterRt , &QPushButton::clicked , mzml_ptr , [=]() mutable {
//       ms1_library_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1Result , this->ui->DoubleSpinBoxMatchMs1MinRt->value() , this->ui->DoubleSpinBoxMatchMs1MaxRt->value());
//    });


//    //一级配对二级
//    connect(this->ui->PushButtonMatchMs1WithMs2 , &QPushButton::clicked , ms_level_matcher_ptr , [=]() mutable{
//        ms_level_matcher_ptr->CopyInfoFromMs1Match(*ms1_library_matcher_ptr);
//        ms_level_matcher_ptr->m_dalton = this->ui->DoubleSpinBoxMatchMs1WithMs2Dalton->value();
//        ms_level_matcher_ptr->m_tolerance_rt_scope = this->ui->DoubleSpinBoxMatchMs1WithMs2ToleranceRt->value();
//        ms_level_matcher_ptr->MatchCardiolipinWithMs2(*mzml_ptr);//配对二级
//        ms_level_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1WithMs2Result);
//        qDebug() << "一级配对二级完成！";
////        //dda
////        ms_level_matcher_ptr->OutPutWithoutMS2WithCsv("D:/qt_project/R_analyze/计算机定性结果/dda");
//    });
//    //根据Rt打印结果
//    connect(this->ui->PushButtonMatchMs1WithMs2FilterRt , &QPushButton::clicked , ms_level_matcher_ptr , [=]() mutable {
//       ms_level_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1WithMs2Result , this->ui->DoubleSpinBoxMatchMs1WithMs2MinRt->value() , this->ui->DoubleSpinBoxMatchMs1WithMs2MaxRt->value());
//    });


//    //找头基
//    connect(this->ui->PushButtonCheakHeadGroup , &QPushButton::clicked , headgroup_finder_ptr , [=]() mutable{
//        headgroup_finder_ptr->CopyInfoFromCardiolipinMatchWithMs2(*ms_level_matcher_ptr);
//        headgroup_finder_ptr->m_ppm = this->ui->DoubleSpinBoxCheakHeadGrouppPPm->value();
//        headgroup_finder_ptr->m_mz_score_weight = this->ui->DoubleSpinBoxCheadHeadGroupMzScoreWeight->value();
//        headgroup_finder_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxCheckHeadGroupPPmWithHalfScore->value();
//        headgroup_finder_ptr->CheckMs2Headgroup();
//        headgroup_finder_ptr->ShowResult(this->ui->TableWidgetCheakHeadGroupResult);
//        qDebug() << "二级找头基完成！！！";
//    });
//    connect(this->ui->PushButtonCheakHeadGroupFilterRt , &QPushButton::clicked , headgroup_finder_ptr , [=]() mutable {
//        headgroup_finder_ptr->ShowResult(this->ui->TableWidgetCheakHeadGroupResult , this->ui->DoubleSpinBoxCheakHeadGroupMinRt->value() , this->ui->DoubleSpinBoxCheakHeadGroupMaxRt->value());
//    });

//    //找PA和FA
//    connect(this->ui->PushButtonFindPaAndFa , &QPushButton::clicked , fragment_finder_ptr , [=]() mutable{
//        fragment_finder_ptr->CopyInfoFromMs2WithHeadgroup(*headgroup_finder_ptr);
//        fragment_finder_ptr->m_ppm = this->ui->DoubleSpinBoxFindPaAndFaPPm->value();
//        fragment_finder_ptr->m_mz_score_weight = this->ui->DoubleSpinBoxFindPaAndFaMzScoreWeight->value();
//        fragment_finder_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxFindPaAndFaPPmWithHalfScore->value();
//        fragment_finder_ptr->FindMs2PaAndFa(*database_ptr);
//        fragment_finder_ptr->ShowResult(this->ui->TableWidgetFindPaAndFaResult);
//        qDebug() << "二级找PA和FA完成！！！！！！";
//    });
//    connect(this->ui->PushButtonFindPaAndFaFilterRt , &QPushButton::clicked , fragment_finder_ptr , [=]() mutable {
//        fragment_finder_ptr->ShowResult(this->ui->TableWidgetFindPaAndFaResult , this->ui->DoubleSpinBoxFindPaAndFaMinRt->value() , this->ui->DoubleSpinBoxFindPaAndFaMaxRt->value());
//    });


//    //拼接
//    connect(this->ui->PushButtonSplice , &QPushButton::clicked , fragment_combiner_ptr , [=]() mutable{
//        fragment_combiner_ptr->CopyInfoFromMs2WithPaAndFa(*fragment_finder_ptr);
//        fragment_combiner_ptr->mode = this->ui->ComboBoxSpliceMode->currentText().toStdString();
//        Cardiolipin::m_delete_redundant_splice_result_radio = 0.2;
//        //增加权重
//        ClSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();
//        MlclSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();
//        DlclSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();

//        ClSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();
//        MlclSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();
//        DlclSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();

//        ClSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();
//        MlclSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();
//        DlclSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();

//        ClSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();
//        MlclSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();
//        DlclSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();
//        Cardiolipin::m_delete_redundant_splice_result_radio = 0.2;

//        fragment_combiner_ptr->splice();
//        fragment_combiner_ptr->ShowSpliceResult(this->ui->TableWidgetSpliceResult);
//        qDebug() << "拼接完成！！！！！！！！！！！！！！！！";

////        //dda
////        fragment_combiner_ptr->OutputResultWithTxt("D:/qt_project/R_analyze/计算机定性结果/dda");//输出txt文件
////        fragment_combiner_ptr->OutPutWithCsv("D:/qt_project/R_analyze/计算机定性结果/dda");//输出csv文件

//    });
//    //splice结果的显示过滤器
//    connect(this->ui->PushButtonSpliceFilterRt , &QPushButton::clicked , fragment_combiner_ptr , [=]( )mutable{
//        fragment_combiner_ptr->ShowSpliceResult(this->ui->TableWidgetSpliceResult , this->ui->DoubleSpinBoxSpliceMinRt->value() , this->ui->DoubleSpinBoxSpliceMaxRt->value());
//    });
//    connect(this->ui->TableWidgetSpliceResult, &QTableWidget::itemClicked, [=](QTableWidgetItem* item){
//        // 获取所选项的行号
//        unsigned int rowIndex = item->row();
//        // 根据行号进行操作
//        fragment_combiner_ptr->PrintSpliceResultByRow(rowIndex);//打印信息
//    });


//    //输出结果
//    connect(this->ui->PushButtonOutputResultDir , &QPushButton::clicked , this , [=](){
//        QFileDialog dialog(this);
//        QString dir_name = dialog.getExistingDirectory();

//        ms_level_matcher_ptr->OutPutWithoutMS2WithCsv(dir_name);//无二级结果

//        fragment_combiner_ptr->OutputResultWithTxt(dir_name);//输出txt文件
//        fragment_combiner_ptr->OutPutWithCsv(dir_name);//输出csv文件

//    });

//    //测试
//    connect(this->ui->PushButtonAddFiles , &QPushButton::clicked , [=](){
//        // dda
//        emit this->SendMs1CsvFileName("D:/qt_project/R_analyze/data/dda/50X_NEG_001(centroid MS1).csv");
//        QStringList file_names = {"D:/qt_project/R_analyze/data/dda/50X_NEG_001(centroid MS2).mzML"};
//         emit this->SendMs2MzmlFileNames(file_names);

//        emit this->ui->PushButtonMatchMs1->clicked();
//        emit this->ui->PushButtonMatchMs1WithMs2->clicked();
//        emit this->ui->PushButtonCheakHeadGroup->clicked();
//        emit this->ui->PushButtonFindPaAndFa->clicked();
//        emit this->ui->PushButtonSplice->clicked();
//    });

//    thread_1->start();

//    emit this->ui->PushButtonAddFiles->clicked();
//}

//idx
MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowTitle(tr("心磷脂定性软件")); // 设置窗口标题为"My Widget"
    thread_1 = new QThread(this);//新建一个线程
    //新建mzml对象，移动到线程thread_1，线程启动
    Mzml* mzml_ptr = new Mzml();
    Database* database_ptr = new Database();

    Ms1LibraryMatcher* ms1_library_matcher_ptr = new Ms1LibraryMatcher();

    MsLevelMatcher* ms_level_matcher_ptr = new MsLevelMatcher();
    HeadgroupFinder* headgroup_finder_ptr = new HeadgroupFinder();
    FragmentFinder* fragment_finder_ptr = new FragmentFinder();
    FragmentCombiner* fragment_combiner_ptr = new FragmentCombiner();
    mzml_ptr->moveToThread(thread_1);
    ms1_library_matcher_ptr->moveToThread(thread_1);
    database_ptr->moveToThread(thread_1);
    ms_level_matcher_ptr->moveToThread(thread_1);
    headgroup_finder_ptr->moveToThread(thread_1);
    fragment_finder_ptr->moveToThread(thread_1);
    fragment_combiner_ptr->moveToThread(thread_1);
   //ui

    this->ui->TableWidgetMatchMs1Result->horizontalHeader()->setStretchLastSection(true);
    this->ui->TableWidgetMatchMs1Result->horizontalHeader()->setVisible(true);
    this->ui->TableWidgetSpliceResult->horizontalHeader()->setVisible(true);

    //打开csv文件，发射信号
    connect(this->ui->ActionOpenMs1 , &QAction::triggered , this , [=](){
        QFileDialog dialog(this);
        dialog.setDirectory("D:/test");
        QString file_name = dialog.getOpenFileName();//获取文件路径
        emit this->SendMs1CsvFileName(file_name);
    });
    //打开多个mzml文件，发射信号
    connect(this->ui->ActionOpenMs2 , &QAction::triggered , this , [=]{
        QFileDialog dialog(this);
        dialog.setDirectory("D:/findmyCL/data/profile");
        QStringList file_names =  dialog.getOpenFileNames();
        emit this->SendMs2MzmlFileNames(file_names);
    });



    //读取一级和二级
    //读取一级csv文件
    connect(this,&MainWindow::SendMs1CsvFileName , mzml_ptr,[=](QString file_name){
        mzml_ptr->ReadMs1FromCsv(file_name);

//        //idx
        mzml_ptr->ConvertMs1RtUnit("mintosec");

        QString message = "Ms1 count:" + QString::number(mzml_ptr->GetLocalMs1Vector().size());
        this->ui->LabelMs1Message->setText(message);

    });
    //读取二级mzml文件，加载数据库
    connect(this,&MainWindow::SendMs2MzmlFileNames , mzml_ptr,[=](QStringList file_names){
        mzml_ptr->ReadMs2FromMzmls(file_names);

//        mzml_ptr->ConvertMs2RtUnit("mintosec");//转换为秒

        database_ptr->LoadAllTable();
        QString message = "Ms2 count:" + QString::number(mzml_ptr->GetLocalMs2Vector().size());
        this->ui->LabelMs2Message->setText(message);

        mzml_ptr->DeleteMs2LowIntensityFragment(0.1);
    });


    //一级配对
    connect(this->ui->PushButtonMatchMs1 , &QPushButton::clicked , mzml_ptr , [=]() mutable {
        ms1_library_matcher_ptr->m_ppm = this->ui->DoubleSpinBoxMatchMs1PPm->value();
        ms1_library_matcher_ptr->m_tolerance_rt_scope = this->ui->DoubleSpinBoxMatchMs1ToleranceRt->value();
        ms1_library_matcher_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxMatchMs1PPmWithHalfScore->value();
        ms1_library_matcher_ptr->MatchMs1WithAllTables(*mzml_ptr , *database_ptr);
        ms1_library_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1Result);//打印结果

    });
    //根据Rt打印结果
    connect(this->ui->PushButtonMatchMs1FilterRt , &QPushButton::clicked , mzml_ptr , [=]() mutable {
       ms1_library_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1Result , this->ui->DoubleSpinBoxMatchMs1MinRt->value() , this->ui->DoubleSpinBoxMatchMs1MaxRt->value());
    });


    //一级配对二级
    connect(this->ui->PushButtonMatchMs1WithMs2 , &QPushButton::clicked , ms_level_matcher_ptr , [=]() mutable{
        ms_level_matcher_ptr->CopyInfoFromMs1Match(*ms1_library_matcher_ptr);
        ms_level_matcher_ptr->m_dalton = this->ui->DoubleSpinBoxMatchMs1WithMs2Dalton->value();
        ms_level_matcher_ptr->m_tolerance_rt_scope = this->ui->DoubleSpinBoxMatchMs1WithMs2ToleranceRt->value();
        ms_level_matcher_ptr->MatchCardiolipinWithMs2(*mzml_ptr);//配对二级
        ms_level_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1WithMs2Result);
        qDebug() << "一级配对二级完成！";
//        //idx
//        ms_level_matcher_ptr->OutPutWithoutMS2WithCsv("D:/qt_project/R_analyze/计算机定性结果/idx");
    });
    //根据Rt打印结果
    connect(this->ui->PushButtonMatchMs1WithMs2FilterRt , &QPushButton::clicked , ms_level_matcher_ptr , [=]() mutable {
       ms_level_matcher_ptr->ShowResult(this->ui->TableWidgetMatchMs1WithMs2Result , this->ui->DoubleSpinBoxMatchMs1WithMs2MinRt->value() , this->ui->DoubleSpinBoxMatchMs1WithMs2MaxRt->value());
    });


    //找头基
    connect(this->ui->PushButtonCheakHeadGroup , &QPushButton::clicked , headgroup_finder_ptr , [=]() mutable{
        headgroup_finder_ptr->CopyInfoFromCardiolipinMatchWithMs2(*ms_level_matcher_ptr);
        headgroup_finder_ptr->m_ppm = this->ui->DoubleSpinBoxCheakHeadGrouppPPm->value();
        headgroup_finder_ptr->m_mz_score_weight = this->ui->DoubleSpinBoxCheadHeadGroupMzScoreWeight->value();
        headgroup_finder_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxCheckHeadGroupPPmWithHalfScore->value();
        headgroup_finder_ptr->CheckMs2Headgroup();
        headgroup_finder_ptr->ShowResult(this->ui->TableWidgetCheakHeadGroupResult);
        qDebug() << "二级找头基完成！！！";
    });
    connect(this->ui->PushButtonCheakHeadGroupFilterRt , &QPushButton::clicked , headgroup_finder_ptr , [=]() mutable {
        headgroup_finder_ptr->ShowResult(this->ui->TableWidgetCheakHeadGroupResult , this->ui->DoubleSpinBoxCheakHeadGroupMinRt->value() , this->ui->DoubleSpinBoxCheakHeadGroupMaxRt->value());
    });

    //找PA和FA
    connect(this->ui->PushButtonFindPaAndFa , &QPushButton::clicked , fragment_finder_ptr , [=]() mutable{
        fragment_finder_ptr->CopyInfoFromMs2WithHeadgroup(*headgroup_finder_ptr);
        fragment_finder_ptr->m_ppm = this->ui->DoubleSpinBoxFindPaAndFaPPm->value();
        fragment_finder_ptr->m_mz_score_weight = this->ui->DoubleSpinBoxFindPaAndFaMzScoreWeight->value();
        fragment_finder_ptr->m_ppm_with_half_score = this->ui->DoubleSpinBoxFindPaAndFaPPmWithHalfScore->value();
        fragment_finder_ptr->FindMs2PaAndFa(*database_ptr);
        fragment_finder_ptr->ShowResult(this->ui->TableWidgetFindPaAndFaResult);
        qDebug() << "二级找PA和FA完成！！！！！！";
    });
    connect(this->ui->PushButtonFindPaAndFaFilterRt , &QPushButton::clicked , fragment_finder_ptr , [=]() mutable {
        fragment_finder_ptr->ShowResult(this->ui->TableWidgetFindPaAndFaResult , this->ui->DoubleSpinBoxFindPaAndFaMinRt->value() , this->ui->DoubleSpinBoxFindPaAndFaMaxRt->value());
    });


    //拼接
    connect(this->ui->PushButtonSplice , &QPushButton::clicked , fragment_combiner_ptr , [=]() mutable{
        fragment_combiner_ptr->CopyInfoFromMs2WithPaAndFa(*fragment_finder_ptr);
        fragment_combiner_ptr->mode = this->ui->ComboBoxSpliceMode->currentText().toStdString();
        Cardiolipin::m_delete_redundant_splice_result_radio = 0.2;
        //增加权重
        ClSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();
        MlclSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();
        DlclSpecificStructure::m_fragment_score_weight = this->ui->DoubleSpinBoxSpliceFragmentScoreWeight->value();

        ClSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();
        MlclSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();
        DlclSpecificStructure::m_pa_exist_score_weight = this->ui->DoubleSpinBoxSplicePaExistScoreWeight->value();

        ClSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();
        MlclSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();
        DlclSpecificStructure::m_fa_consistency_score_weight = this->ui->DoubleSpinBoxSpliceFaConsistencyScoreWeight->value();


        ClSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();
        MlclSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();
        DlclSpecificStructure::m_fa_intensity_variance_score_weight = this->ui->DoubleSpinBoxSpliceFaIntensityVarianceScoreWeight->value();


        fragment_combiner_ptr->splice();
        fragment_combiner_ptr->ShowSpliceResult(this->ui->TableWidgetSpliceResult);
        qDebug() << "拼接完成！！！！！！！！！！！！！！！！";
    });
    //splice结果的显示过滤器
    connect(this->ui->PushButtonSpliceFilterRt , &QPushButton::clicked , fragment_combiner_ptr , [=]( )mutable{
        fragment_combiner_ptr->ShowSpliceResult(this->ui->TableWidgetSpliceResult , this->ui->DoubleSpinBoxSpliceMinRt->value() , this->ui->DoubleSpinBoxSpliceMaxRt->value());
    });
    connect(this->ui->TableWidgetSpliceResult, &QTableWidget::itemClicked, [=](QTableWidgetItem* item){
        // 获取所选项的行号
        unsigned int rowIndex = item->row();
        // 根据行号进行操作
        fragment_combiner_ptr->PrintSpliceResultByRow(rowIndex);//打印信息
    });

    //输出结果
    connect(this->ui->PushButtonOutputResultDir , &QPushButton::clicked , this , [=](){
        QFileDialog dialog(this);
        QString dir_name = dialog.getExistingDirectory();

        ms_level_matcher_ptr->OutPutWithoutMS2WithCsv(dir_name);//无二级结果

        fragment_combiner_ptr->OutputResultWithTxt(dir_name);//输出txt文件
        fragment_combiner_ptr->OutPutWithCsv(dir_name);//输出csv文件
    });

    //测试
    connect(this->ui->PushButtonAddFiles , &QPushButton::clicked , [=](){
        // IDX
        emit this->SendMs1CsvFileName("D:/test/peaks _copy.csv");
        QStringList file_names = {"D:/findmyCL/data/profile/NEG_ID_01.mzML(centroid)","D:/findmyCL/data/profile/NEG_ID_02.mzML(centroid)","D:/findmyCL/data/profile/NEG_ID_03.mzML(centroid)","D:/findmyCL/data/profile/NEG_ID_04.mzML(centroid)"};
        //QStringList file_names = {"D:/findmyCL/data/profile/NEG_ID_01.mzML","D:/findmyCL/data/profile/NEG_ID_02.mzML","D:/findmyCL/data/profile/NEG_ID_03.mzML","D:/findmyCL/data/profile/NEG_ID_04.mzML"};
       // QStringList file_names = {"D:/findmyCL/data/profile/NEG_ID_01.mzML(centroid)"};
        emit this->SendMs2MzmlFileNames(file_names);
        emit this->ui->PushButtonMatchMs1->clicked();
        emit this->ui->PushButtonMatchMs1WithMs2->clicked();
        emit this->ui->PushButtonCheakHeadGroup->clicked();
        emit this->ui->PushButtonFindPaAndFa->clicked();
        emit this->ui->PushButtonSplice->clicked();
    });

    thread_1->start();

    emit this->ui->PushButtonAddFiles->clicked();
}

MainWindow::~MainWindow()
{
    delete ui;
    this->thread_1->quit();
}

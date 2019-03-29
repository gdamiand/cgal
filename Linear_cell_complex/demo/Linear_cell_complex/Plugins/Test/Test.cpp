#include <QTime>
#include <QApplication>
#include <QAction>
#include <QMainWindow>
#include <QMessageBox>
#include <QStringList>

#include "MainWindow.h"
#include "LCC_demo_plugin_interface.h"

using namespace CGAL;

class LCC_demo_test_plugin :
  public QObject,
  public LCC_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::LCC_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "LCC_Demo.PluginInterface/1.0")
  
public:
  void init(MainWindow *mymw, Scene *myscene)
  {
    scene = myscene;
    mw = mymw;
    
    QAction *myAction = new QAction("Test plugin", mw);
    myAction->setProperty("subMenuName","Test me");
    connect(myAction, SIGNAL(triggered()), this, SLOT(on_actionTest_triggered()));
    _actions <<myAction;
  }
  
  QList<QAction*> actions()const {return _actions;}

public Q_SLOTS:
  void on_actionTest_triggered();
  
private:
  Scene* scene;
  MainWindow* mw;
  QList<QAction*> _actions;
};

void LCC_demo_test_plugin::on_actionTest_triggered()
{
  QMessageBox msgBox;
  msgBox.setText("Test worked.");
  msgBox.exec();
  Q_EMIT(mw->sceneChanged());
}

#include "Test.moc"

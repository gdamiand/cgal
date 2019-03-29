#ifndef LCC_DEMO_PLUGIN_INTERFACE_H
#define LCC_DEMO_PLUGIN_INTERFACE_H

#include <QString>
#include <QList>
#include <QtPlugin>
#include <QDebug>
#include "typedefs.h"

class QAction;
class QMainWindow;
class MainWindow;
namespace CGAL
{
    /*!
     * This virtual class provides the basic functions used for making a plugin.
     */
    class LCC_demo_plugin_interface
    {
    public:
        //! \brief Initializes the plugin
        //! This function acts like a constructor. This is where the attributes must be initialized.
        //! The Message_interface allows to print warnings or errors on the screen and the `Console` widget.
        virtual void init(MainWindow*, Scene*) = 0;

        //!Contains all the plugin's actions.
        virtual QList<QAction*> actions() const = 0;
    };
}

Q_DECLARE_INTERFACE(CGAL::LCC_demo_plugin_interface, "LCC_Demo.PluginInterface/1.0")

#endif // LCC_DEMO_PLUGIN_INTERFACE_H

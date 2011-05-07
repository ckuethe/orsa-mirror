TEMPLATE = subdirs

unix:!macx { SUBDIRS += orsaPDS.dynamiclib.pro }
SUBDIRS += orsaPDS.staticlib.pro

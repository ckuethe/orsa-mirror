#include <orsaPDS/RadioScienceGravity.h>

#include <iostream>
#include <cstdio>

using namespace orsaPDS;

RadioScienceGravityFile::RadioScienceGravityFile(const std::string & fileName) {
    data = new RadioScienceGravityData;
    FILE * fp = fopen(fileName.c_str(),"r");
    if (!fp) return;
    fclose(fp);
}

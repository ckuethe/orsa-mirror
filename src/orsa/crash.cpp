#include <orsa/crash.h>

#include <orsa/debug.h>

#include <stdlib.h>

// for alternatives to orsa::crash() see i.e. http://code.google.com/p/google-coredumper/

void orsa::crash() {
    ORSA_ERROR("time to die!");
    abort();
}

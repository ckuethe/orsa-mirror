#include <orsa/debug.h>
#include <cstdlib>
#include <libgen.h>
#include <cstdio>

int main(int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    double D    = 20.0; // diameter, km
    double step =  2.5;
    
    const double R=0.5*D;

    const bool inside_crater_only = true;
    
    // before all, write this rule
    // with this, all intermediate files are kept instead of deleted automatically
    printf(".SECONDARY:\n");
    printf("\n");

    {
        // first, the targets
        printf("all:");
        double x=-R;
        while (x<=R) {
            double y=-R;
            while (y<=R) {
                
                if ((inside_crater_only) && (x*x+y*y<=R*R)) {
                    // printf(" TSC_history_%+.3f_%+.3f.out",x,y);
                    printf(" TSC_resume_%+.3f_%+.3f.out",x,y);
                }
                
                y += step;
            }
            x += step;
        }
    }
    
    printf("\n"); 
    printf("\n");
    
    {
        // now the rules
        double x=-R;
        while (x<=R) {
            double y=-R;
            while (y<=R) {
                
                if ((inside_crater_only) && (x*x+y*y<=R*R)) {
                    printf("TSC_resume_%+.3f_%+.3f.out:\n",x,y,x,y);
                    printf("\t./ThermalStressCraters %+.3f %+.3f > $*.log 2>&1\n",x,y);
                    printf("\n");
                }
                
                y += step;
            }
            x += step;
        }
    }
    
    return 0;
}

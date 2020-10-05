#include "simplex.h"
#include "SH2ijk.h"

#include <orsa/legendre.h>
#include <orsa/util.h>
#include <orsa/statistic.h>

#include "shape.h"

#include <gsl/gsl_sf_legendre.h>

// norm_coeff = normalization_factor * coeff
double normalization_factor(const size_t & l,
                            const size_t & m) {
    return orsa::normalization_sphericalHarmonicsToNormalizedSphericalHarmonics(l,m).get_d();
}

double radius(const double & theta,
              const double & phi,
              const SHcoeff & norm_A,
              const SHcoeff & norm_B) {
    const double c_theta = cos(theta);
    double r=0;
    if (norm_A.size() != norm_B.size()) {
        ORSA_DEBUG("problems...");
    }
    const double max_degree = norm_A.size()-1;
    double * Plm_array = new double[gsl_sf_legendre_array_n(max_degree)];
    gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,max_degree,c_theta,+1,Plm_array);
    for (size_t l=0; l<norm_A.size(); ++l) {
        if (norm_A[l].size() != norm_B[l].size()) {
            ORSA_DEBUG("problems...");
        }
        for (size_t m=0; m<norm_A[l].size(); ++m) {
            if ( (norm_A[l][m] != 0.0) ||
                 (norm_B[l][m] != 0.0) ) {
                const double Nf = normalization_factor(l,m);
                // const double Plm = orsa::LegendreP(l,m,c_theta).get_d();
                const double Plm = Plm_array[gsl_sf_legendre_array_index(l,m)];
                // ORSA_DEBUG("LegendreP[%i][%i](%g) = %g",l,m,c_theta,Plm);
                r += Plm*(norm_A[l][m]*cos(m*phi))/Nf;
                if (m!=0) {
                    r += Plm*(norm_B[l][m]*sin(m*phi))/Nf;
                }   
            }
        }
    }
    delete[] Plm_array;
    return r;
}

double radius(orsa::Vector & r, const SHcoeff & norm_A, const SHcoeff & norm_B) {
    const orsa::Vector u = r.normalized();
    const double theta = acos(u.getZ());
    const double phi   = atan2(u.getY(),u.getX());
    //
    r = u*radius(theta,phi,norm_A,norm_B);
    return r.length();
}

int main(int argc, char **argv) {
    
    /*
#warning comment out!! testing only...
    ORSA_DEBUG("WARNING: using fixed random seed in production!!!");
    // test specific cases, for debug purposes only!
    orsa::GlobalRNG::randomSeed = 7777;
    */
    
    const double km = orsa::FromUnits(1,orsa::Unit::KM);
    
    if (argc != 6) {
        printf("Usage: %s <plate-model-file-in> <SH-file-out> <max-degree> <epsabs-km> <epsrel>\n",argv[0]);
        exit(0);
    }
    
    const std::string  inputFile = argv[1];
    const std::string outputFile = argv[2];
    const size_t max_degree = strtoul(argv[3],0,10);
    const double epsabs     = orsa::FromUnits(fabs(atof(argv[4])),orsa::Unit::KM);
    const double epsrel     = fabs(atof(argv[5]));
    
    
    osg::ref_ptr<InputShape> shapeModel = new InputShape;
       if (!shapeModel->read(inputFile)) {
       ORSA_ERROR("problems encountered while reading shape file...");
       exit(0);
    }
        
    /*
    osg::ref_ptr<GaskellPlateModel> shapeModel = new GaskellPlateModel;
    if (!shapeModel->read(inputFile)) {
        ORSA_ERROR("problems encountered while reading shape file...");
        exit(0);
    }
    */
    
    // init normalized coefficients
    SHcoeff norm_A;
    SHcoeff norm_B;
    norm_A.resize(max_degree+1);
    norm_B.resize(max_degree+1);
    for (size_t l=0; l<=max_degree; ++l) {
        norm_A[l].resize(l+1);
        norm_B[l].resize(l+1);
        for (size_t m=0; m<=l; ++m) {
            norm_A[l][m] = 0;
            norm_B[l][m] = 0;
        }
    }
    
    // estimated coefficients
    std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > stat_norm_A, stat_norm_B;
    stat_norm_A.resize(max_degree+1);
    stat_norm_B.resize(max_degree+1);
    for (size_t l=0; l<=max_degree; ++l) {
        stat_norm_A[l].resize(l+1);
        stat_norm_B[l].resize(l+1);
        for (size_t m=0; m<=l; ++m) {
            stat_norm_A[l][m] = new orsa::Statistic<double>;
            stat_norm_B[l][m] = new orsa::Statistic<double>;
        }
    }
    
    double * Plm_array = new double[gsl_sf_legendre_array_n(max_degree)];
    
    // iter
    size_t counter=0;
    while (1) {
        ++counter;
        
        double dx,dy,dz;
        GlobalRNG::instance()->rng()->gsl_ran_dir_3d(&dx,&dy,&dz);
        const double theta = acos(dz);
        const double phi   = atan2(dy,dx);
        //
        orsa::Vector shape_intersectionPoint;
        orsa::Vector shape_normal;
        shapeModel->rayIntersection(shape_intersectionPoint,
                                    shape_normal,
                                    orsa::Vector(0,0,0),
                                    orsa::Vector(dx,dy,dz),
                                    false);
        const double r = shape_intersectionPoint.length();
        
        const double c_theta = cos(theta);
        
        gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE,max_degree,c_theta,+1,Plm_array);
        
        for (size_t l=0; l<=max_degree; ++l) {
            for (size_t m=0; m<=l; ++m) {
                // const double Plm = LegendreP(l,m,c_theta).get_d();
                const double Plm = Plm_array[gsl_sf_legendre_array_index(l,m)];
                stat_norm_A[l][m]->insert(r*Plm*cos(m*phi));
                if (m!=0) stat_norm_B[l][m]->insert(r*Plm*sin(m*phi));
                // ORSA_DEBUG("TEST: %g == %g ??",Plm,LegendreP(l,m,c_theta).get_d());
            }
        }
        
        {
            static size_t jj=0;
            ++jj;
            if (jj%10000==0) {
                bool done=true;
                for (size_t l=0; l<=max_degree; ++l) {
                    for (size_t m=0; m<=l; ++m) {
                        const double Nf = normalization_factor(l,m);
                        const double factor = mpf_class(mpq_class((2*l+1)*(2-orsa::kronecker(0,m))*orsa::factorial(l-m),(orsa::factorial(l+m)))).get_d();
                        // ORSA_DEBUG("factor[%i][%i] = %16.6f",l,m,factor);
                        
                        // save
                        norm_A[l][m] = factor*Nf*stat_norm_A[l][m]->average();
                        if (m!=0) norm_B[l][m] = factor*Nf*stat_norm_B[l][m]->average();
                        // output file
                        writeSH(norm_A,norm_B,outputFile,orsa::Unit::KM);
                        
                        double val, sigma;
                        char msg[4096];
                        //
                        val   = factor*Nf*stat_norm_A[l][m]->average();
                        sigma = factor*Nf*stat_norm_A[l][m]->averageError();
                        //
                        bool bad=false;
                        if (sigma>epsabs) bad=true;
                        if ((fabs(val)*epsrel<sigma) && (fabs(val)>3*sigma)) {
                            bad=true;
                        }
                        if (bad) {
                            done=false;
                            sprintf(msg,"(not converged)");
                        } else {
                            sprintf(msg,"");
                        }
                        //
                        ORSA_DEBUG("stat_norm_A[%02i][%02i]->average() = %16.6f +/- %16.6f   %s",l,m,
                                   factor*Nf*stat_norm_A[l][m]->average()/km,
                                   factor*Nf*stat_norm_A[l][m]->averageError()/km,
                                   msg);
                        
                        if (m!=0) {
                            //
                            val   = factor*Nf*stat_norm_B[l][m]->average();
                            sigma = factor*Nf*stat_norm_B[l][m]->averageError();
                            //
                            bool bad=false;
                            if (sigma>epsabs) bad=true;
                            if ((fabs(val)*epsrel<sigma) && (fabs(val)>3*sigma)) {
                                bad=true;
                            }
                            if (bad) {
                                done=false;
                                sprintf(msg,"(not converged)");
                            } else {
                                sprintf(msg,"");
                            }
                            //
                            ORSA_DEBUG("stat_norm_B[%02i][%02i]->average() = %16.6f +/- %16.6f   %s",l,m,
                                             factor*Nf*stat_norm_B[l][m]->average()/km,
                                             factor*Nf*stat_norm_B[l][m]->averageError()/km,
                                             msg);
                        }
                    }
                }
                if (done) break;
            }
        }
    }
    
    delete[] Plm_array;
    
    return 0;
}

#ifndef _ORSA_INPUT_OUTPUT_MPC_
#define _ORSA_INPUT_OUTPUT_MPC_

#include <orsa/datetime.h>

#include <string>

namespace orsaInputOutput {
  
    int         MPC_charToInt(const char c);
    std::string MPC_intToChar(const int i);
  
    orsa::Time  MPC_packedToTime(const std::string & packedEpoch);
    std::string MPC_timeToPacked(const orsa::Time & t);
  
    unsigned int MPC_packedNumber(const std::string & packedNumber);
    std::string  MPC_packNumber(const unsigned int & number, const size_t & digits);
    
    std::string MPC_packedDesignation(const std::string & readableDesignation);
    
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_

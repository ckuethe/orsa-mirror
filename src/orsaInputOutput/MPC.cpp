#include <orsaInputOutput/MPC.h>

using namespace orsa;
using namespace orsaInputOutput;

#include <orsaSolarSystem/datetime.h>

int orsaInputOutput::MPC_charToInt(const char c) {
    int d;
    const int ch = (int)c;
    switch (ch) {
        case '0': d=0;  break;
        case '1': d=1;  break;
        case '2': d=2;  break;
        case '3': d=3;  break;
        case '4': d=4;  break;
        case '5': d=5;  break;
        case '6': d=6;  break;
        case '7': d=7;  break;
        case '8': d=8;  break;
        case '9': d=9;  break;
        case 'A': d=10; break; 
        case 'B': d=11; break;
        case 'C': d=12; break;
        case 'D': d=13; break;
        case 'E': d=14; break;
        case 'F': d=15; break;
        case 'G': d=16; break;
        case 'H': d=17; break;
        case 'I': d=18; break;
        case 'J': d=19; break;
        case 'K': d=20; break;
        case 'L': d=21; break;
        case 'M': d=22; break;
        case 'N': d=23; break;
        case 'O': d=24; break;
        case 'P': d=25; break;
        case 'Q': d=26; break;
        case 'R': d=27; break;
        case 'S': d=28; break;
        case 'T': d=29; break;
        case 'U': d=30; break;
        case 'V': d=31; break;
        case 'W': d=32; break;
        case 'X': d=33; break;
        case 'Y': d=34; break;
        case 'Z': d=35; break;
        case 'a': d=36; break;
        case 'b': d=37; break;
        case 'c': d=38; break;
        case 'd': d=39; break;
        case 'e': d=40; break;
        case 'f': d=41; break;
        case 'g': d=42; break;
        case 'h': d=43; break;
        case 'i': d=44; break;
        case 'j': d=45; break;
        case 'k': d=46; break;
        case 'l': d=47; break;
        case 'm': d=48; break;
        case 'n': d=49; break;
        case 'o': d=50; break;
        case 'p': d=51; break;
        case 'q': d=52; break;
        case 'r': d=53; break;
        case 's': d=54; break;
        case 't': d=55; break;
        case 'u': d=56; break;
        case 'v': d=57; break;
        case 'w': d=58; break;
        case 'x': d=59; break;
        case 'y': d=60; break;
        case 'z': d=61; break;               
        default:
            ORSA_DEBUG("case not handled, c: [%c]",c);
            d=0;
    }
    return d;
}

std::string orsaInputOutput::MPC_intToChar(const int i) {
    char c;
    switch (i) {
        case  0: c='0'; break;
        case  1: c='1'; break;
        case  2: c='2'; break;
        case  3: c='3'; break;
        case  4: c='4'; break;
        case  5: c='5'; break;
        case  6: c='6'; break;
        case  7: c='7'; break;
        case  8: c='8'; break;
        case  9: c='9'; break;
        case 10: c='A'; break;
        case 11: c='B'; break;
        case 12: c='C'; break;
        case 13: c='D'; break;
        case 14: c='E'; break;
        case 15: c='F'; break;
        case 16: c='G'; break;
        case 17: c='H'; break;
        case 18: c='I'; break;
        case 19: c='J'; break;
        case 20: c='K'; break;
        case 21: c='L'; break;
        case 22: c='M'; break;
        case 23: c='N'; break;
        case 24: c='O'; break;
        case 25: c='P'; break;
        case 26: c='Q'; break;
        case 27: c='R'; break;
        case 28: c='S'; break;
        case 29: c='T'; break;
        case 30: c='U'; break;
        case 31: c='V'; break;
        case 32: c='W'; break;
        case 33: c='X'; break;
        case 34: c='Y'; break;
        case 35: c='Z'; break;
        case 36: c='a'; break;
        case 37: c='b'; break;
        case 38: c='c'; break;
        case 39: c='d'; break;
        case 40: c='e'; break;
        case 41: c='f'; break;
        case 42: c='g'; break;
        case 43: c='h'; break;
        case 44: c='i'; break;
        case 45: c='j'; break;
        case 46: c='k'; break;
        case 47: c='l'; break;
        case 48: c='m'; break;
        case 49: c='n'; break;
        case 50: c='o'; break;
        case 51: c='p'; break;
        case 52: c='q'; break;
        case 53: c='r'; break;
        case 54: c='s'; break;
        case 55: c='t'; break;
        case 56: c='u'; break;
        case 57: c='v'; break;
        case 58: c='w'; break;
        case 59: c='x'; break;
        case 60: c='y'; break;
        case 61: c='z'; break;
        default:
            ORSA_DEBUG("case not handled, i: [%i]",i);
            c=' ';
            orsa::crash();
    }
    char str[1024];
    sprintf(str,"%c",c);
    return str;
}

orsa::Time orsaInputOutput::MPC_packedToTime(const std::string & packedEpoch) {
    const char * s = packedEpoch.c_str();
    if (strlen(s) != 5) {
        ORSA_DEBUG("problems");
        return orsa::Time();
    }
    int y, m, d;
    y  = 100*MPC_charToInt(s[0]);
    y +=  10*MPC_charToInt(s[1]);
    y +=     MPC_charToInt(s[2]);
    m  = MPC_charToInt(s[3]);
    d  = MPC_charToInt(s[4]);
    const orsa::Time t = 
        orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d,0,0,0,0), 
                                       orsaSolarSystem::TS_TDT);
    return t;
}

std::string orsaInputOutput::MPC_timeToPacked(const orsa::Time & t_in) {
    const orsa::Time t_TDT = orsaSolarSystem::ToTimeScale(t_in,orsaSolarSystem::TS_TDT);
    int y,m,d,H,M,S,ms;
    orsaSolarSystem::gregorDay(t_TDT,y,m,d,H,M,S,ms);
    if (H || M || S || ms) {
        ORSA_DEBUG("problems...");
    }
    char pkd[1024];
    pkd[0] = '\0';
    strcat(pkd,MPC_intToChar(y/100).c_str());
    strcat(pkd,MPC_intToChar((y%100)/10).c_str());
    strcat(pkd,MPC_intToChar(y%10).c_str());
    strcat(pkd,MPC_intToChar(m).c_str());
    strcat(pkd,MPC_intToChar(d).c_str());
    return pkd;
}

unsigned int orsaInputOutput::MPC_packedNumber(const std::string & packedNumber) {
    if (strlen(packedNumber.c_str()) == 0) {
        return 0;
    }
    if ((strlen(packedNumber.c_str()) != 2) && (strlen(packedNumber.c_str()) != 5)) {
        // ORSA_DEBUG("MPC packed number use only 5 digits (arg: [%s])",packedNumber.c_str());
        // but compact designations use 2 digits compact numbers too...
        return 0;
    }
    // if last digit is a 'P', issue an error, as it is probably a comet number
    if (packedNumber[strlen(packedNumber.c_str())-1] == 'P') {
        // ORSA_DEBUG("[%s] interpreted as comet number, returning zero.",packedNumber.c_str());
        return 0;
    }
    unsigned int num=0;
    for (unsigned int k=0; k<strlen(packedNumber.c_str()); ++k) {
        if (!isalnum(packedNumber[k])) {
            ORSA_DEBUG("cannot handle this character: [%s]",packedNumber[k]);
        }
        num *= 10;
        num += MPC_charToInt(packedNumber[k]);
    }
    // ORSA_DEBUG("[%s] -> %i",packedNumber.c_str(),num);
    return num;
}

std::string orsaInputOutput::MPC_packNumber(const unsigned int & number, const size_t & digits) {
    unsigned int p10=1;
    for (size_t k=0; k<digits; ++k) {
        p10 *= 10;
    }
    unsigned int n10=1;
    unsigned int d10=0;
    while (n10<number) {
        n10 *= 10;
        d10 += 1;
    }
    if (number<p10) {
        // trivially return the number
        char fmt[4096];
        switch (digits) {
            case 0: sprintf(fmt,"%%00i"); break;
            case 1: sprintf(fmt,"%%01i"); break;
            case 2: sprintf(fmt,"%%02i"); break;
            case 3: sprintf(fmt,"%%03i"); break;
            case 4: sprintf(fmt,"%%04i"); break;
            case 5: sprintf(fmt,"%%05i"); break;
            case 6: sprintf(fmt,"%%06i"); break;
            case 7: sprintf(fmt,"%%07i"); break;
            case 8: sprintf(fmt,"%%08i"); break;
            case 9: sprintf(fmt,"%%09i"); break;
            default:
            ORSA_DEBUG("problems...");
            exit(0);
        }
        char pkd[4096];
        sprintf(pkd,fmt,number);
        return pkd;
    }
    std::vector<std::string> digit;
    digit.resize(digits);
    char full[4096];
    sprintf(full,"%u",number);
    // ORSA_DEBUG("number: %u   digits: %u   p10: %i   n10: %i   d10: %i   full: [%s]",number,digits,p10,n10,d10,full);
    for (size_t k=0; k<digits; ++k) {
        if (k!=(digits-1)) {
            digit[k] = full[d10-k-1];
        } else {
            digit[k] = orsaInputOutput::MPC_intToChar(number/(p10/10));
        }
        // ORSA_DEBUG("digit[%u]: [%s]",k,digit[k].c_str());
    }
    char pkd[4096];
    sprintf(pkd,"");
    for (size_t k=0; k<digits; ++k) {
        strcat(pkd,digit[digits-k-1].c_str());
    }
    return pkd;
}

std::string orsaInputOutput::MPC_packedDesignation(const std::string & readableDesignation) {
    // input is a string like "2013 GJ69" and output is "K13G69J"
    // there are exceptions such as "2040 P-L", see http://www.minorplanetcenter.net/iau/info/PackedDes.html
    char p1[4096], p2[4096]; // two parts
    sscanf(readableDesignation.c_str(),"%s %s",p1,p2);
    if ((strlen(p1)==0) || (strlen(p2)==0)) {
        return readableDesignation;
    }
    const std::string s1(p1);
    const std::string s2(p2);
    if (s2=="P-L") {
        char pkd[4096];
        sprintf(pkd,"PLS%s",p1);
        return pkd;
    } else if (s2=="T-1") {
        char pkd[4096];
        sprintf(pkd,"T1S%s",p1);
        return pkd;
    } else if (s2=="T-2") {
        char pkd[4096];
        sprintf(pkd,"T2S%s",p1);
        return pkd;
    } else if (s2=="T-3") {
        char pkd[4096];
        sprintf(pkd,"T3S%s",p1);
        return pkd;
    } else {
        // standard
        const int year = atoi(p1);
        if ((year<1800) || (year>2100)) {
            // probable composite name, i.e. "El Djezair"
            // ORSA_DEBUG("problem with [%s]",readableDesignation.c_str());
            return readableDesignation;
        }
        char c1, c2;
        int num=0;
        sscanf(p2,"%c%c%i",&c1,&c2,&num);
        char pkd[4096];
        sprintf(pkd,"%s%02i%c%s%c",orsaInputOutput::MPC_intToChar(year/100).c_str(),year%100,c1,orsaInputOutput::MPC_packNumber(num,2).c_str(),c2);
        return pkd;
    }
}


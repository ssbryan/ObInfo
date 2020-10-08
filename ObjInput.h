
// ObjInput.h : header file
//

#pragma once

#include <string>
#include <vector>
//#include <cmath>

struct Header
{
    Header(void)
    {
        time = 0;
        nbodies = 0;
        ndim = 0;
        nsph = 0;
        ndark = 0;
        nstar = 0;
    }
    
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
    int pad;
};

struct GasParticle
{
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float eps;
    float metals ;
    float phi;
};

struct DarkParticle
{
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
};

struct StarParticle
{
    float mass;
    float pos[3];
    float vel[3];
    float metals ;
    float tform;
    float eps;
    float phi;
};

struct GasParticleD
{
    double mass;
    double pos[3];
    double vel[3];
    double rho;
    double temp;
    double eps;
    double metals ;
    double phi;
};

struct DarkParticleD
{
    double mass;
    double pos[3];
    double vel[3];
    double eps;
    double phi;
};

struct StarParticleD
{
    double mass;
    double pos[3];
    double vel[3];
    double metals ;
    double tform;
    double eps;
    double phi;
};

class ObjInput
{
public:
    ObjInput(std::string& inf, std::string& outf);
    bool    ReadAndWriteInfo();
    bool    ReadObjects(void);
    void    WriteInfo(void);

private:
    std::vector<GasParticle> gparts;
    std::vector<DarkParticle> dparts;
    std::vector<StarParticle> sparts;
    Header  h;

    std::string mInFName;
    std::string mOutFName;
};

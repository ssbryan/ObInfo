
// ObGenDlg.cpp : implementation file
//

#include "ObjInput.h"

#include <fstream>
#include <assert.h>


ObjInput::ObjInput(std::string& inf, std::string& outf)
: mInFName(inf)
, mOutFName(outf)
{
}

bool ObjInput::ReadAndWriteInfo(void)
{	
    bool okay = ReadObjects();

    if (okay)
    {
        WriteInfo();
    }
    else
    {
        printf("failed to read file\n");
        return false;
    }

    return true;
}

bool ObjInput::ReadObjects(void)
{
    std::ifstream infile;
    infile.open(mInFName, std::ios::in | std::ios::binary);
    bool flipped = false;

    if (infile.is_open())
    {
        bool okay = true;
        
        // binary format: header, nsph gas particles, ndark dark particles, nstar star particles

        infile.read((char*)(&h), sizeof(Header));

        if (h.ndim != 3)
        {
            unsigned char flipbuf[4];
            unsigned char* ch = (unsigned char*)&(h.ndim);

            for (int i = 0; i < 4; ++i)
            {
                flipbuf[3 - i] = ch[i];
            }

            int* pdim = (int*)flipbuf;;
            int dim = *pdim;

            if (dim == 3)
            {
                flipped = true;
                h.ndim = dim;

                // switch the other vals and fix the header
                int vals[4] = {h.nbodies, h.nsph, h.ndark, h.nstar};
                int flipvals[4];

                for (int j = 0; j < 4; ++j)
                {
                    ch = (unsigned char*)&(vals[j]);

                    for (int i = 0; i < 4; ++i)
                    {
                        flipbuf[3 - i] = ch[i];
                    }

                    int* pnum = (int*)flipbuf;
                    flipvals[j] = *pnum;
                }

                h.nbodies = flipvals[0];
                h.nsph = flipvals[1];
                h.ndark = flipvals[2];
                h.nstar = flipvals[3];
            }
        }
    }

    int start = infile.tellg();
    infile.seekg (0, infile.end);
    int fileLen = infile.tellg();
    int dataLen = fileLen - start;
    infile.seekg(start);

    int tryFloat = h.nsph * sizeof(GasParticle) + 
        h.ndark * sizeof(DarkParticle) +
        h.nstar * sizeof(StarParticle);

    int tryDouble = h.nsph * sizeof(GasParticleD) + 
        h.ndark * sizeof(DarkParticleD) +
        h.nstar * sizeof(StarParticleD);

    bool useF = (dataLen == tryFloat);
    bool useD = (dataLen == tryDouble);
    assert(useF && !useD);
    
    for (int i = 0; i < h.nsph; ++i)
    {
        GasParticle g;
        char gbuf[sizeof(GasParticle)];

        if (!flipped)
        {
            infile.read((char*)(&g), sizeof(GasParticle));
        }
        else
        {
            // flip each 4 byte chunk
            infile.read((char*)(gbuf), sizeof(GasParticle));
            char swap[sizeof(GasParticle)];

            for (int j = 0; j < sizeof(GasParticle); j += 4)
            {
                swap[j + 3] = gbuf[j + 0];
                swap[j + 2] = gbuf[j + 1];
                swap[j + 1] = gbuf[j + 2];
                swap[j + 0] = gbuf[j + 3];
            }

            float* nums = (float*)swap;
            g.mass = nums[0];
            g.pos[0] = nums[1];
            g.pos[1] = nums[2];
            g.pos[2] = nums[3];
            g.vel[0] = nums[4];
            g.vel[1] = nums[5];
            g.vel[2] = nums[6];
            g.rho = nums[7];
            g.temp = nums[8];
            g.eps = nums[9];
            g.metals = nums[10];
            g.phi = nums[11];
        }
        
        gparts.push_back(g);
    }

    for (int i = 0; i < h.ndark; ++i)
    {
        DarkParticle d;
        char dbuf[sizeof(DarkParticle)];
        
        if (!flipped)
        {
            infile.read((char*)(&d), sizeof(DarkParticle));
        }
        else
        {
            // flip each 4 byte chunk
            infile.read((char*)(dbuf), sizeof(DarkParticle));
            char swap[sizeof(DarkParticle)];

            for (int j = 0; j < sizeof(DarkParticle); j += 4)
            {
                swap[j + 3] = dbuf[j + 0];
                swap[j + 2] = dbuf[j + 1];
                swap[j + 1] = dbuf[j + 2];
                swap[j + 0] = dbuf[j + 3];
            }

            float* nums = (float*)swap;
            d.mass = nums[0];
            d.pos[0] = nums[1];
            d.pos[1] = nums[2];
            d.pos[2] = nums[3];
            d.vel[0] = nums[4];
            d.vel[1] = nums[5];
            d.vel[2] = nums[6];
            d.eps = nums[7];
            d.phi = nums[8];
        }
        
        dparts.push_back(d);
    }

    for (int i = 0; i < h.nstar; ++i)
    {
        StarParticle s;
        char sbuf[sizeof(StarParticle)];
    
        if (!flipped)
        {
            infile.read((char*)(&s), sizeof(StarParticle));
        }
        else
        {
            // flip each 4 byte chunk
            infile.read((char*)(sbuf), sizeof(StarParticle));
            char swap[sizeof(StarParticle)];

            for (int j = 0; j < sizeof(StarParticle); j += 4)
            {
                swap[j + 3] = sbuf[j + 0];
                swap[j + 2] = sbuf[j + 1];
                swap[j + 1] = sbuf[j + 2];
                swap[j + 0] = sbuf[j + 3];
            }

            float* nums = (float*)swap;
            s.mass = nums[0];
            s.pos[0] = nums[1];
            s.pos[1] = nums[2];
            s.pos[2] = nums[3];
            s.vel[0] = nums[4];
            s.vel[1] = nums[5];
            s.vel[2] = nums[6];
            s.metals = nums[7];
            s.tform = nums[8];
            s.eps = nums[9];
            s.phi = nums[10];
        }
            
        sparts.push_back(s);
    }

    infile.close();
    return true;
}

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157E+308
#endif

void ObjInput::WriteInfo(void)
{
    // mass range
    double minMass = DBL_MAX;
    double maxMass = 0;

    // loc extents
    double minX = DBL_MAX;
    double minY = DBL_MAX;
    double minZ = DBL_MAX;
    double maxX = -DBL_MAX;
    double maxY = -DBL_MAX;
    double maxZ = -DBL_MAX;

    // vel extents
    double minVX = DBL_MAX;
    double minVY = DBL_MAX;
    double minVZ = DBL_MAX;
    double maxVX = -DBL_MAX;
    double maxVY = -DBL_MAX;
    double maxVZ = -DBL_MAX;

    // get statistics for loc, vel
    double avgx = 0.0;
    double avgy = 0.0;
    double avgz = 0.0;
    double avgvx = 0.0;
    double avgvy = 0.0;
    double avgvz = 0.0;

    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;
    double totmass = 0.0;
    double totmomx = 0.0;
    double totmomy = 0.0;
    double totmomz = 0.0;

    // gas particle properties
    double rhoMaxG = 0;
    double rhoMinG = 1e20;
    double rhoAvgG = 0;
    double tempMaxG = 0;
    double tempMinG = 1e20;
    double tempAvgG = 0;
    double epsMaxG = 0;
    double epsMinG = 1e20;
    double epsAvgG = 0;
    double metalsMaxG = 0;
    double metalsMinG = 1e20;
    double metalsAvgG = 0;
    double phiMaxG = 0;
    double phiMinG = 1e20;
    double phiAvgG = 0;

    // dark particle properties
    double epsMaxD = 0;
    double epsMinD = 1e20;
    double epsAvgD = 0;
    double phiMaxD = 0;
    double phiMinD = 1e20;
    double phiAvgD = 0;

    // star particle properties
    double metalsMaxS = 0;
    double metalsMinS = 1e20;
    double metalsAvgS = 0;
    double tformMaxS = 0;
    double tformMinS = 1e20;
    double tformAvgS = 0;
    double epsMaxS = 0;
    double epsMinS = 1e20;
    double epsAvgS = 0;
    double phiMaxS = 0;
    double phiMinS = 1e20;
    double phiAvgS = 0;

    // traverse particle vectors and create string

    if (h.nsph > 0)
    {
        int gnum = gparts.size();
        
        for (std::vector<GasParticle>::iterator gIt = gparts.begin();
             gIt != gparts.end();
             ++gIt)
        {
            double x = gIt->pos[0]; 
            double y = gIt->pos[1]; 
            double z = gIt->pos[2]; 
            double vx = gIt->vel[0]; 
            double vy = gIt->vel[1]; 
            double vz = gIt->vel[2];

            if (x < minX)
            {
                minX = x;
            }
            
            if (x > maxX)
            {
                maxX = x;
            }
            
            if (y < minY)
            {
                minY = y;
            }
            
            if (y > maxY)
            {
                maxY = y;
            }
            
            if (z < minZ)
            {
                minZ = z;
            }
            
            if (z > maxZ)
            {
                maxZ = z;
            }
            
            if (vx < minVX)
            {
                minVX = vx;
            }
            
            if (vx > maxVX)
            {
                maxVX = vx;
            }
            
            if (vy < minVY)
            {
                minVY = vy;
            }
            
            if (vy > maxVY)
            {
                maxVY = vy;
            }
            
            if (vz < minVZ)
            {
                minVZ = vz;
            }
            
            if (vz > maxVZ)
            {
                maxVZ = vz;
            }
            
	        avgx += x;
	        avgy += y;
	        avgz += z;
	        avgvx += vx;
	        avgvy += vy;
	        avgvz += vz;

            // get CoM and total momentum
	        double m = gIt->mass;
	        xc += m * x;
	        yc += m * y;
	        zc += m * z;

            if (m < minMass)
            {
                minMass = m;
            }

            if (m > maxMass)
            {
                maxMass = m;
            }
            
	        totmass += m;

	        totmomx += m * vx;
	        totmomy += m * vy;
	        totmomz += m * vz;

            // properties
            if (gIt->rho > rhoMaxG)
            {
                rhoMaxG = gIt->rho;
            }

            if (gIt->rho < rhoMinG)
            {
                rhoMinG = gIt->rho;
            }

            rhoAvgG += gIt->rho;

            if (gIt->temp > tempMaxG)
            {
                tempMaxG = gIt->temp;
            }

            if (gIt->temp < tempMinG)
            {
                tempMinG = gIt->temp;
            }

            tempAvgG += gIt->temp;

            if (gIt->eps > epsMaxG)
            {
                epsMaxG = gIt->eps;
            }

            if (gIt->eps < epsMinG)
            {
                epsMinG = gIt->eps;
            }

            epsAvgG += gIt->eps;

            if (gIt->metals > metalsMaxG)
            {
                metalsMaxG = gIt->metals;
            }

            if (gIt->metals < metalsMinG)
            {
                metalsMinG = gIt->metals;
            }

            metalsAvgG += gIt->metals;

            if (gIt->phi > phiMaxG)
            {
                phiMaxG = gIt->phi;
            }

            if (gIt->phi < phiMinG)
            {
                phiMinG = gIt->phi;
            }

            phiAvgG += gIt->phi;
        }

        rhoAvgG /= h.nsph;
        tempAvgG /= h.nsph;
        epsAvgG /= h.nsph;
        metalsAvgG /= h.nsph;
        phiAvgG /= h.nsph;
    }
    
    if (h.ndark > 0)
    {
        for (std::vector<DarkParticle>::iterator dIt = dparts.begin();
             dIt != dparts.end();
             ++dIt)
        {
            double x = dIt->pos[0]; 
            double y = dIt->pos[1]; 
            double z = dIt->pos[2]; 
            double vx = dIt->vel[0]; 
            double vy = dIt->vel[1]; 
            double vz = dIt->vel[2];
            
            if (x < minX)
            {
                minX = x;
            }
            
            if (x > maxX)
            {
                maxX = x;
            }
            
            if (y < minY)
            {
                minY = y;
            }
            
            if (y > maxY)
            {
                maxY = y;
            }
            
            if (z < minZ)
            {
                minZ = z;
            }
            
            if (z > maxZ)
            {
                maxZ = z;
            }
            
            if (vx < minVX)
            {
                minVX = vx;
            }
            
            if (vx > maxVX)
            {
                maxVX = vx;
            }
            
            if (vy < minVY)
            {
                minVY = vy;
            }
            
            if (vy > maxVY)
            {
                maxVY = vy;
            }
            
            if (vz < minVZ)
            {
                minVZ = vz;
            }
            
            if (vz > maxVZ)
            {
                maxVZ = vz;
            }
            
	        avgx += x;
	        avgy += y;
	        avgz += z;
	        avgvx += vx;
	        avgvy += vy;
	        avgvz += vz;

            // get CoM and total momentum
	        double m = dIt->mass;
	        xc += m * x;
	        yc += m * y;
	        zc += m * z;

            if (m < minMass)
            {
                minMass = m;
            }

            if (m > maxMass)
            {
                maxMass = m;
            }
            
	        totmass += m;

	        totmomx += m * vx;
	        totmomy += m * vy;
	        totmomz += m * vz;

            //properties
            if (dIt->eps > epsMaxD)
            {
                epsMaxD = dIt->eps;
            }

            if (dIt->eps < epsMinD)
            {
                epsMinD = dIt->eps;
            }

            epsAvgD += dIt->eps;

            if (dIt->phi > phiMaxD)
            {
                phiMaxD = dIt->phi;
            }

            if (dIt->phi < phiMinD)
            {
                phiMinD = dIt->phi;
            }

            phiAvgD += dIt->phi;
        }

        epsAvgD /= h.ndark;
        phiAvgD /= h.ndark;
    }
    
    if (h.nstar > 0)
    {
        for (std::vector<StarParticle>::iterator sIt = sparts.begin();
             sIt != sparts.end();
             ++sIt)
        {
            double x = sIt->pos[0]; 
            double y = sIt->pos[1]; 
            double z = sIt->pos[2]; 
            double vx = sIt->vel[0]; 
            double vy = sIt->vel[1]; 
            double vz = sIt->vel[2];
            
            if (x < minX)
            {
                minX = x;
            }
            
            if (x > maxX)
            {
                maxX = x;
            }
            
            if (y < minY)
            {
                minY = y;
            }
            
            if (y > maxY)
            {
                maxY = y;
            }
            
            if (z < minZ)
            {
                minZ = z;
            }
            
            if (z > maxZ)
            {
                maxZ = z;
            }
            
            if (vx < minVX)
            {
                minVX = vx;
            }
            
            if (vx > maxVX)
            {
                maxVX = vx;
            }
            
            if (vy < minVY)
            {
                minVY = vy;
            }
            
            if (vy > maxVY)
            {
                maxVY = vy;
            }
            
            if (vz < minVZ)
            {
                minVZ = vz;
            }
            
            if (vz > maxVZ)
            {
                maxVZ = vz;
            }
            
	        avgx += x;
	        avgy += y;
	        avgz += z;
	        avgvx += vx;
	        avgvy += vy;
	        avgvz += vz;

            // get CoM and total momentum
	        double m = sIt->mass;
	        xc += m * x;
	        yc += m * y;
	        zc += m * z;

            if (m < minMass)
            {
                minMass = m;
            }

            if (m > maxMass)
            {
                maxMass = m;
            }
            
	        totmass += m;

	        totmomx += m * vx;
	        totmomy += m * vy;
	        totmomz += m * vz;

            // properties
            if (sIt->tform > tformMaxS)
            {
                tformMaxS = sIt->tform;
            }

            if (sIt->tform < tformMinS)
            {
                tformMinS = sIt->tform;
            }

            tformAvgS += sIt->tform;

            if (sIt->eps > epsMaxS)
            {
                epsMaxS = sIt->eps;
            }

            if (sIt->eps < epsMinS)
            {
                epsMinS = sIt->eps;
            }

            epsAvgS += sIt->eps;

            if (sIt->metals > metalsMaxS)
            {
                metalsMaxS = sIt->metals;
            }

            if (sIt->metals < metalsMinS)
            {
                metalsMinS = sIt->metals;
            }

            metalsAvgS += sIt->metals;

            if (sIt->phi > phiMaxS)
            {
                phiMaxS = sIt->phi;
            }

            if (sIt->phi < phiMinS)
            {
                phiMinS = sIt->phi;
            }

            phiAvgS += sIt->phi;
        }

        tformAvgS /= h.nstar;
        epsAvgS /= h.nstar;
        metalsAvgS /= h.nstar;
        phiAvgS /= h.nstar;
    }
    

    avgx /= h.nbodies;
    avgy /= h.nbodies;
    avgz /= h.nbodies;
    avgvx /= h.nbodies;
    avgvy /= h.nbodies;
    avgvz /= h.nbodies;

    // analyze input data to get
    // --total momentum
    // --total angular momentum
    // --CoM compared to 0,0,0
    // --total kinetic energy

    xc /= totmass;
    yc /= totmass;
    zc /= totmass;

    // use cross product: m * cross(v, r)
    // to get angular momentum
    double angmomx = 0.0;
    double angmomy = 0.0;
    double angmomz = 0.0;

    if (h.nsph > 0)
    {
        for (std::vector<GasParticle>::iterator gIt = gparts.begin();
             gIt != gparts.end();
             ++gIt)
        {
            double rx = gIt->pos[0] - xc;
            double ry = gIt->pos[1] - yc;
            double rz = gIt->pos[2] - zc;
            double crossx = gIt->vel[1] * rz - gIt->vel[2] * ry;
            double crossy = gIt->vel[2] * rx - gIt->vel[0] * rz;
            double crossz = gIt->vel[0] * ry - gIt->vel[1] * rx;
            angmomx += gIt->mass * crossx;
            angmomy += gIt->mass * crossy;
            angmomz += gIt->mass * crossz;
        }
    }

    if (h.ndark > 0)
    {
        for (std::vector<DarkParticle>::iterator dIt = dparts.begin();
             dIt != dparts.end();
             ++dIt)
        {
            double rx = dIt->pos[0] - xc;
            double ry = dIt->pos[1] - yc;
            double rz = dIt->pos[2] - zc;
            double crossx = dIt->vel[1] * rz - dIt->vel[2] * ry;
            double crossy = dIt->vel[2] * rx - dIt->vel[0] * rz;
            double crossz = dIt->vel[0] * ry - dIt->vel[1] * rx;
            angmomx += dIt->mass * crossx;
            angmomy += dIt->mass * crossy;
            angmomz += dIt->mass * crossz;
        }
    }

    if (h.nstar > 0)
    {
        for (std::vector<StarParticle>::iterator sIt = sparts.begin();
             sIt != sparts.end();
             ++sIt)
        {
            double rx = sIt->pos[0] - xc;
            double ry = sIt->pos[1] - yc;
            double rz = sIt->pos[2] - zc;
            double crossx = sIt->vel[1] * rz - sIt->vel[2] * ry;
            double crossy = sIt->vel[2] * rx - sIt->vel[0] * rz;
            double crossz = sIt->vel[0] * ry - sIt->vel[1] * rx;
            angmomx += sIt->mass * crossx;
            angmomy += sIt->mass * crossy;
            angmomz += sIt->mass * crossz;
        }
    }

    std::ofstream ofsa;

    ofsa.open(mOutFName, std::ios::out);
 
    ofsa << "Number of objects: " << h.nbodies << std::endl;
    ofsa << "Number of SPH objects: " << h.nsph << std::endl;
    ofsa << "Number of Dark objects: " << h.ndark << std::endl;
    ofsa << "Number of Star objects: " << h.nstar << std::endl;
    ofsa << "Mass: min " << minMass << ", max " << maxMass << ", total " << totmass << ", avg " << totmass / h.nbodies << std::endl;
    ofsa << "Center of Mass: " << xc << ", " << yc << ", " << zc << std::endl;
    ofsa << "Average mass: " << totmass / h.nbodies << std::endl;
    ofsa << "Location range: X " << minX << " to " << maxX << ", Y " << minY << " to " << maxY << ", Z " << minZ << " to " << maxZ << std::endl;
    ofsa << "Average location: " << avgx << ", " << avgy << ", " << avgz << std::endl;
    ofsa << "Velocity range: X " << minVX << " to " << maxVX << ", Y " << minVY << " to " << maxVY << ", Z " << minVZ << " to " << maxVZ << std::endl;
    ofsa << "Average velocity: " << avgvx << ", " << avgvy << ", " << avgvz << std::endl;
    ofsa << "Momentum: " << totmomx << ", " << totmomy << ", " << totmomz << std::endl;

    if (h.nsph > 0)
    {
        ofsa << "SPH Properties: " << std::endl;
        ofsa << "    rho: min " << rhoMinG << ", max " << rhoMaxG << ", avg " << rhoAvgG << std::endl;
        ofsa << "    temp: min " << tempMinG << ", max " << tempMaxG << ", avg " << tempAvgG << std::endl;
        ofsa << "    eps: min " << epsMinG << ", max " << epsMaxG << ", avg " << epsAvgG << std::endl;
        ofsa << "    metals: min " << metalsMinG << ", max " << metalsMaxG << ", avg " << metalsAvgG << std::endl;
        ofsa << "    phi: min " << phiMinG << ", max " << phiMaxG << ", avg " << phiAvgG << std::endl;
    }

    if (h.ndark > 0)
    {
        ofsa << "Dark Properties: " << std::endl;
        ofsa << "    eps: min " << epsMinD << ", max " << epsMaxD << ", avg " << epsAvgD << std::endl;
        ofsa << "    phi: min " << phiMinD << ", max " << phiMaxD << ", avg " << phiAvgD << std::endl;
    }

    if (h.nstar > 0)
    {
        ofsa << "Star Properties: " << std::endl;
        ofsa << "    tform: min " << tformMinS << ", max " << tformMaxS << ", avg " << tformAvgS << std::endl;
        ofsa << "    eps: min " << epsMinS << ", max " << epsMaxS << ", avg " << epsAvgS << std::endl;
        ofsa << "    metals: min " << metalsMinS << ", max " << metalsMaxS << ", avg " << metalsAvgS << std::endl;
        ofsa << "    phi: min " << phiMinS << ", max " << phiMaxS << ", avg " << phiAvgS << std::endl;
    }

    ofsa << "Ang Momentum: " << angmomx << ", " << angmomy << ", " << angmomz << std::endl;

    ofsa.close();
}
